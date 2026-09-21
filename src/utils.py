import pandas as pd
import pyranges as pr

def standardize_chrom_strand(df):
    """Standardizes chromosome prefix and strand format for rMATS tables."""
    df['Chromosome'] = df['chr']
    if not df['Chromosome'].str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    df['Strand'] = df['strand']
    return df

def merge_overlap_genomic(lc_df, rmats_introns_df, max_boundary_drift=30):
    """Performs strand-aware PyRanges spatial join with a maximum boundary drift threshold."""
    pr_lc = pr.PyRanges(lc_df)
    pr_rmats = pr.PyRanges(rmats_introns_df)
    
    overlap_df = pr_rmats.join(pr_lc, suffix='_LeafCutter').df
    
    if overlap_df.empty:
        return overlap_df
        
    overlap_df['Start_Drift'] = (overlap_df['Start'] - overlap_df['Start_LeafCutter']).abs()
    overlap_df['End_Drift'] = (overlap_df['End'] - overlap_df['End_LeafCutter']).abs()
    
    filtered_matches = overlap_df[
        (overlap_df['Start_Drift'] <= max_boundary_drift) & 
        (overlap_df['End_Drift'] <= max_boundary_drift)
    ].copy()
    
    return filtered_matches