import pandas as pd
import pyranges as pr
from scipy.stats import mannwhitneyu

JUNCTION_KEY = ['Chromosome', 'Strand', 'Start', 'End']
TOOLS = ['LeafCutter', 'rMATS', 'MAJIQ']

# Column names of each tool's junction coordinates in the merged table
MERGED_COLUMNS = {
    'LeafCutter': ['Chromosome', 'Strand', 'Start', 'End'],
    'rMATS': ['Chromosome', 'Strand_rMATS', 'Start_rMATS', 'End_rMATS'],
    'MAJIQ': ['Chromosome', 'Strand_MAJIQ', 'Start_MAJIQ', 'End_MAJIQ'],
}

IN_MERGED = 'In merged data'
MATCHES_OTHER = 'Not merged: matches another tool'
NO_MATCH = 'Not merged: no match in other tools'
STATUSES = [IN_MERGED, MATCHES_OTHER, NO_MATCH]


def unique_junctions(df):
    """One row per distinct junction (rMATS repeats a junction across events)."""
    return df[JUNCTION_KEY].drop_duplicates().reset_index(drop=True)


def _with_partner(a, b, max_boundary_drift):
    """Boolean Series over a's rows: True if b has a junction within the drift on the same strand."""
    left = a.rename_axis('idx').reset_index()
    joined = pr.PyRanges(left).join(pr.PyRanges(b), suffix='_b').df
    if joined.empty:
        return pd.Series(False, index=a.index)
    close = joined[
        ((joined['Start'] - joined['Start_b']).abs() <= max_boundary_drift)
        & ((joined['End'] - joined['End_b']).abs() <= max_boundary_drift)
    ]
    return pd.Series(a.index.isin(close['idx']), index=a.index)


def _unmerged_detail(tool, junctions, merged, max_boundary_drift):
    """Unique junctions of one tool with a flag for presence in the merged data and, for the ones
    missing from it, whether each other tool has a counterpart (columns 'has_<tool>')."""
    mine = junctions[tool]
    in_merged = merged[MERGED_COLUMNS[tool]].drop_duplicates()
    in_merged.columns = JUNCTION_KEY
    flagged = mine.merge(in_merged, on=JUNCTION_KEY, how='left', indicator=True)
    mine = mine.assign(merged=(flagged['_merge'] == 'both').values)

    missing = mine[~mine['merged']]
    for other in (t for t in TOOLS if t != tool):
        mine[f'has_{other}'] = _with_partner(missing, junctions[other], max_boundary_drift).reindex(
            mine.index, fill_value=False)
    return mine


def classify_junctions(tool, tool_data, merged, max_boundary_drift):
    """Unique junctions of `tool` with a 'Status' column: merged, unmerged with a counterpart in
    another tool, or unmerged with no counterpart anywhere."""
    junctions = {t: unique_junctions(df) for t, df in tool_data.items()}
    detail = _unmerged_detail(tool, junctions, merged, max_boundary_drift)
    others = [f'has_{t}' for t in TOOLS if t != tool]
    has_any = detail[others].any(axis=1)
    detail['Status'] = NO_MATCH
    detail.loc[has_any, 'Status'] = MATCHES_OTHER
    detail.loc[detail['merged'], 'Status'] = IN_MERGED
    return detail[JUNCTION_KEY + ['Status']]


def unmerged_junctions(tool_data, merged, max_boundary_drift):
    """Counts, per tool, the filtered junctions that are absent from the final merged data.

    Each unmerged junction is put in exactly one category: no counterpart in either other tool,
    a counterpart in only one of them, or counterparts in both (so it was dropped by the
    closest-match step rather than for lack of a partner). Returns one row per tool.
    """
    junctions = {t: unique_junctions(df) for t, df in tool_data.items()}

    rows = []
    for tool in TOOLS:
        detail = _unmerged_detail(tool, junctions, merged, max_boundary_drift)
        missing = detail[~detail['merged']]
        others = [t for t in TOOLS if t != tool]
        a, b = (missing[f'has_{o}'] for o in others)

        row = {
            'Tool': tool,
            'Filtered': len(detail),
            'In_Merged': int(detail['merged'].sum()),
            'Not_In_Merged': len(missing),
            'No_Match_In_Either_Other': int((~a & ~b).sum()),
        }
        for t in TOOLS:  # 0 for the tool itself, so every row has the same columns
            row[f'Only_Matches_{t}'] = 0
        row[f'Only_Matches_{others[0]}'] = int((a & ~b).sum())
        row[f'Only_Matches_{others[1]}'] = int((b & ~a).sum())
        row['Matches_Both_But_Dropped'] = int((a & b).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def rmats_rate_comparison(rmats, tool_data, merged, max_boundary_drift):
    """Compares rMATS read rates of junctions that did / did not reach the merged data.

    Rates are per event: mean over replicates of IJC / IncFormLen (inclusion rate) and
    SJC / SkipFormLen (skipping rate), averaged over the two groups. Per unique junction the event
    with the highest rate is used. Events whose count lists are malformed (not one value per
    replicate) have no rate and are left out, and the number left out is reported.

    Returns (summary DataFrame with one row per status, {'Inc': {...}, 'Skip': {...}} of
    Mann-Whitney p-values for each unmerged status versus the merged junctions).
    """
    status = classify_junctions('rMATS', tool_data, merged, max_boundary_drift)

    events = rmats.assign(
        Inc_Rate=(rmats['Inc_Rate_Hyp'] + rmats['Inc_Rate_Norm']) / 2,
        Skip_Rate=(rmats['Skip_Rate_Hyp'] + rmats['Skip_Rate_Norm']) / 2,
    )
    per_junction = events.groupby(JUNCTION_KEY, as_index=False)[['Inc_Rate', 'Skip_Rate']].max()
    data = per_junction.merge(status, on=JUNCTION_KEY)  # max() skips NaN, so all-NaN stays NaN

    rows = []
    pvalues = {'Inc': {}, 'Skip': {}}
    reference = data[(data['Status'] == IN_MERGED) & data['Inc_Rate'].notna()]
    for st in STATUSES:
        sub_all = data[data['Status'] == st]
        sub = sub_all[sub_all['Inc_Rate'].notna()]
        row = {'Status': st, 'Junctions': len(sub_all), 'Excluded_Malformed': len(sub_all) - len(sub)}
        for kind in ('Inc', 'Skip'):
            col = f'{kind}_Rate'
            row[f'Median_{kind}_Rate'] = sub[col].median()
            row[f'Q1_{kind}_Rate'] = sub[col].quantile(0.25)
            row[f'Q3_{kind}_Rate'] = sub[col].quantile(0.75)
            if st != IN_MERGED and len(sub) and len(reference):
                pvalues[kind][st] = mannwhitneyu(sub[col].dropna(), reference[col].dropna(),
                                                 alternative='two-sided').pvalue
        rows.append(row)
    return pd.DataFrame(rows), pvalues
