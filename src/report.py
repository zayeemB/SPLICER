import numpy as np

TOOLS = ['LeafCutter', 'rMATS', 'MAJIQ']


def _pct(n, total):
    return f"{n} ({100 * n / total:.1f}%)" if total else "0"


def sign_agreement(summary):
    """Returns a boolean Series: True where all three tools report the same DeltaPSI sign."""
    signs = np.sign(summary[[f'{t}_DeltaPSI' for t in TOOLS]])
    return signs.nunique(axis=1) == 1


def write_report(summary, filter_counts, report_path, disagreements_path, max_boundary_drift,
                 unmerged=None, rate_summary=None, rate_pvalues=None):
    """Writes a short plain-text sign-agreement report and a CSV of the disagreeing rows."""
    delta = {t: summary[f'{t}_DeltaPSI'] for t in TOOLS}
    signs = {t: np.sign(s) for t, s in delta.items()}
    agree = sign_agreement(summary)
    total = len(summary)

    lines = ["SPLICER multi-tool report", "========================", ""]

    lines.append("Rows passing per-tool filters (before merging)")
    for tool, (n_in, n_out, rule) in filter_counts.items():
        lines.append(f"  {tool:<11} {n_out:>6} / {n_in:<7} {rule}")
    lines += [f"  Matching: strand-aware, within {max_boundary_drift} bp, closest match only",
              f"  Merged consensus rows: {total}", ""]

    if unmerged is not None:
        lines.append("Filtered junctions not in the final merged data (unique junctions)")
        for _, r in unmerged.iterrows():
            lines.append(f"  {r['Tool']:<11} {r['Not_In_Merged']:>5} of {r['Filtered']:<5} "
                         f"({100 * r['Not_In_Merged'] / r['Filtered']:.1f}%) not merged")
            lines.append(f"      no match in either other tool: {r['No_Match_In_Either_Other']}")
            for t in TOOLS:
                if t != r['Tool']:
                    lines.append(f"      only matches {t}: {r[f'Only_Matches_{t}']}")
            lines.append(f"      matches both but dropped by closest-match step: {r['Matches_Both_But_Dropped']}")
        lines.append("")

    if rate_summary is not None:
        lines.append("rMATS read rates of junctions in vs. not in the merged data")
        lines.append("  Rate = mean over replicates of reads / effective length, averaged over both groups;")
        lines.append("  best event per junction. Events with malformed count lists are excluded.")
        for kind, title in (('Inc', 'Inclusion rate (IJC / IncFormLen)'), ('Skip', 'Skipping rate (SJC / SkipFormLen)')):
            lines.append(f"  {title}")
            lines.append(f"    {'':<38}{'junctions':>9}  {'excluded':>8}  {'median (IQR)':>20}")
            for _, r in rate_summary.iterrows():
                iqr = f"{r[f'Median_{kind}_Rate']:.3f} ({r[f'Q1_{kind}_Rate']:.3f}-{r[f'Q3_{kind}_Rate']:.3f})"
                lines.append(f"    {r['Status']:<38}{r['Junctions']:>9}  {r['Excluded_Malformed']:>8}  {iqr:>20}")
            for status, p in (rate_pvalues or {}).get(kind, {}).items():
                lines.append(f"    Mann-Whitney vs merged, '{status}': p = {p:.2g}")
        lines.append("")

    lines.append("DeltaPSI sign agreement (Hyp - Norm)")
    lines.append(f"  All three agree:  {_pct(int(agree.sum()), total)}")
    lines.append(f"  Any disagreement: {_pct(int((~agree).sum()), total)}")
    for i, a in enumerate(TOOLS):
        for b in TOOLS[i + 1:]:
            same = int((signs[a] == signs[b]).sum())
            lines.append(f"  {a} vs {b}: {_pct(same, total)}")
    lines.append("")

    lines.append("Odd one out (disagreements where only one tool differs)")
    for t in TOOLS:
        a, b = [u for u in TOOLS if u != t]
        odd = (signs[a] == signs[b]) & (signs[t] != signs[a])
        lines.append(f"  {t:<11} {int(odd.sum())}")
    lines.append("")

    lines.append("By event type (rows / sign disagreements)")
    for ev, idx in summary.groupby('EventType').groups.items():
        lines.append(f"  {ev:<5} {len(idx):>5} / {int((~agree.loc[idx]).sum())}")
    lines.append("")

    lines.append("Spearman correlation of DeltaPSI")
    for i, a in enumerate(TOOLS):
        for b in TOOLS[i + 1:]:
            lines.append(f"  {a} vs {b}: {delta[a].corr(delta[b], method='spearman'):.2f}")
    lines += ["", f"Rows with disagreeing signs are listed in {disagreements_path.name}"]

    report_path.write_text("\n".join(lines) + "\n")
    summary[~agree].to_csv(disagreements_path, index=False)
