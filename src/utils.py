import pandas as pd
import pyranges as pr


def standardize_chrom_strand(df):
    """Standardizes chromosome prefix and strand format for rMATS tables."""
    df['Chromosome'] = df['chr']
    if not df['Chromosome'].str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    df['Strand'] = df['strand']
    return df


def merge_overlap_genomic(df1, df2, max_boundary_drift=2, suffix='_right', closest_only=True):
    """
    Performs strand-aware PyRanges spatial join with a maximum boundary drift threshold,
    using a dynamic suffix to prevent column collisions during sequential joins.

    If closest_only is True, each row of df1 keeps only its closest match(es) in df2 (smallest
    summed start+end drift; exact ties are all kept), so a junction is not paired with a
    neighbouring junction that merely falls within the tolerance.
    """
    df1 = df1.copy()
    df1['_row'] = range(len(df1))
    pr_1 = pr.PyRanges(df1)
    pr_2 = pr.PyRanges(df2)

    # Perform spatial join using the dynamic suffix
    overlap_df = pr_1.join(pr_2, suffix=suffix).df

    if overlap_df.empty:
        return overlap_df.drop(columns='_row', errors='ignore')

    # Boundary drift per match, named by suffix so sequential merges don't overwrite each other
    start_drift, end_drift = f'Start_Drift{suffix}', f'End_Drift{suffix}'
    overlap_df[start_drift] = (overlap_df['Start'] - overlap_df[f'Start{suffix}']).abs()
    overlap_df[end_drift] = (overlap_df['End'] - overlap_df[f'End{suffix}']).abs()

    filtered_matches = overlap_df[
        (overlap_df[start_drift] <= max_boundary_drift) &
        (overlap_df[end_drift] <= max_boundary_drift)
    ].copy()

    if closest_only and not filtered_matches.empty:
        total_drift = filtered_matches[start_drift] + filtered_matches[end_drift]
        best = total_drift.groupby(filtered_matches['_row']).transform('min')
        filtered_matches = filtered_matches[total_drift == best]

    return filtered_matches.drop(columns='_row')
