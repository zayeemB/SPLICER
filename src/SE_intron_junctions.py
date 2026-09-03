import pandas as pd
import pyranges as pr

def parse_leafcutter(leafcutter_path):
    """
    Parses LeafCutter intron effect sizes/significance table.
    Expected 'intron' column format: chr:start:end:clu_X_strand or chr:start:end:clu_X
    """
    lc_df = pd.read_csv(leafcutter_path, sep=r'\s+')
    
    split_coords = lc_df['intron'].str.split(':', expand=True)
    lc_df['Chromosome'] = split_coords[0]
    lc_df['Start'] = split_coords[1].astype(int)
    
    # Subtract 1 to convert LeafCutter's 1-based end to 0-based half-open coordinate
    lc_df['End'] = split_coords[2].astype(int) - 1
    
    return lc_df

def parse_rmats_se(rmats_path):
    """
    Parses rMATS SE (Skipped Exon) output and generates matching intron intervals:
    - Upstream intron: upstreamEE to exonStart_0base
    - Downstream intron: exonEnd to downstreamES
    - Combined intron (exon exclusion): upstreamEE to downstreamES
    """
    rmats_df = pd.read_csv(rmats_path, sep='\t')
    
    # Standardize chromosome column name
    rmats_df['Chromosome'] = rmats_df['chr']
    if not rmats_df['Chromosome'].str.startswith('chr').any():
        rmats_df['Chromosome'] = 'chr' + rmats_df['Chromosome'].astype(str)
        
    # Generate long-format table for all introns associated with each SE event
    intron_rows = []
    for idx, row in rmats_df.iterrows():
        event_id = row['ID']
        gene_symbol = row['geneSymbol']
        strand = row['strand']
        chr_name = row['Chromosome']
        
        # Upstream intron
        intron_rows.append({
            'rMATS_ID': event_id,
            'geneSymbol': gene_symbol,
            'Chromosome': chr_name,
            'Start': row['upstreamEE'],
            'End': row['exonStart_0base'],
            'Intron_Type': 'Upstream_Intron',
            'rMATS_IncLevelDifference': row['IncLevelDifference'],
            'rMATS_FDR': row['FDR']
        })
        # Downstream intron
        intron_rows.append({
            'rMATS_ID': event_id,
            'geneSymbol': gene_symbol,
            'Chromosome': chr_name,
            'Start': row['exonEnd'],
            'End': row['downstreamES'],
            'Intron_Type': 'Downstream_Intron',
            'rMATS_IncLevelDifference': row['IncLevelDifference'],
            'rMATS_FDR': row['FDR']
        })
        # Skipped/Exclusion intron
        intron_rows.append({
            'rMATS_ID': event_id,
            'geneSymbol': gene_symbol,
            'Chromosome': chr_name,
            'Start': row['upstreamEE'],
            'End': row['downstreamES'],
            'Intron_Type': 'Exclusion_Intron',
            'rMATS_IncLevelDifference': row['IncLevelDifference'],
            'rMATS_FDR': row['FDR']
        })

    return pd.DataFrame(intron_rows)

def merge_exact_introns(lc_df, rmats_introns_df):
    """
    Merges datasets based on exact intron boundary coordinates (Chromosome, Start, End).
    """
    merged = pd.merge(
        rmats_introns_df,
        lc_df,
        on=['Chromosome', 'Start', 'End'],
        how='inner',
        suffixes=('_rMATS', '_LeafCutter')
    )
    return merged

def merge_overlap_genomic(lc_df, rmats_introns_df):
    """
    Merges datasets using PyRanges to capture overlapping intervals 
    even if exact splice site boundaries differ slightly.
    """
    pr_lc = pr.PyRanges(lc_df)
    pr_rmats = pr.PyRanges(rmats_introns_df)
    
    # Join overlapping ranges
    overlap = pr_rmats.join(pr_lc, suffix='_LeafCutter').df
    return overlap

# Example Execution
if __name__ == "__main__":
    # File paths
    leafcutter_file = "leafcutter_effect_sizes.txt"
    rmats_se_file = "SE.MATS.JC.txt"
    
    # Load and format inputs
    lc_data = parse_leafcutter(leafcutter_file)
    rmats_introns = parse_rmats_se(rmats_se_file)
    
    # Strategy 1: Exact Intron Boundary Matching
    exact_matches = merge_exact_introns(lc_data, rmats_introns)
    exact_matches.to_csv("rmats_leafcutter_exact_matches.tsv", sep="\t", index=False)
    
    # Strategy 2: Spatial Overlap Matching
    overlap_matches = merge_overlap_genomic(lc_data, rmats_introns)
    overlap_matches.to_csv("rmats_leafcutter_overlap_matches.tsv", sep="\t", index=False)