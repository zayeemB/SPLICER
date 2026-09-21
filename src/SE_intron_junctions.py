import pandas as pd
import pyranges as pr

def parse_leafcutter(leafcutter_path):
    """
    Parses LeafCutter CSV output.
    Extracts Chromosome, Start (0-based), End (converted to 0-based), and Strand.
    """
    lc_df = pd.read_csv(leafcutter_path)
    
    # Split intron string (Format: chr:start:end:clu_X_strand or chr:start:end:clu_X)
    split_coords = lc_df['intron'].str.split(':', expand=True)
    lc_df['Chromosome'] = split_coords[0]
    lc_df['Start'] = split_coords[1].astype(int)
    
    # Subtract 1 to convert LeafCutter 1-based end to 0-based half-open coordinate
    lc_df['End'] = split_coords[2].astype(int) - 1
    
    # Extract strand from cluster identifier (e.g., clu_100_- -> '-')
    lc_df['Strand'] = split_coords[3].str.split('_').str[-1]
    
    # Standardize chromosome formatting (ensure 'chr' prefix)
    if not lc_df['Chromosome'].str.startswith('chr').any():
        lc_df['Chromosome'] = 'chr' + lc_df['Chromosome'].astype(str)
        
    return lc_df

def parse_rmats_se(rmats_path):
    """
    Parses rMATS SE CSV output into a long-format intron table.
    """
    rmats_df = pd.read_csv(rmats_path)
    
    # Standardize chromosome and strand column names for PyRanges
    rmats_df['Chromosome'] = rmats_df['chr']
    if not rmats_df['Chromosome'].str.startswith('chr').any():
        rmats_df['Chromosome'] = 'chr' + rmats_df['Chromosome'].astype(str)
        
    rmats_df['Strand'] = rmats_df['strand']
        
    intron_rows = []
    for idx, row in rmats_df.iterrows():
        event_id = row['ID']
        gene_symbol = row.get('geneSymbol', '')
        strand_val = row['Strand']
        chr_name = row['Chromosome']
        
        # Upstream intron
        intron_rows.append({
            'rMATS_ID': event_id,
            'geneSymbol': gene_symbol,
            'Chromosome': chr_name,
            'Strand': strand_val,
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
            'Strand': strand_val,
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
            'Strand': strand_val,
            'Start': row['upstreamEE'],
            'End': row['downstreamES'],
            'Intron_Type': 'Exclusion_Intron',
            'rMATS_IncLevelDifference': row['IncLevelDifference'],
            'rMATS_FDR': row['FDR']
        })

    return pd.DataFrame(intron_rows)

def merge_exact_introns(lc_df, rmats_introns_df):
    """
    1:1 Exact Coordinate & Strand Match.
    """
    merged = pd.merge(
        rmats_introns_df,
        lc_df,
        on=['Chromosome', 'Start', 'End', 'Strand'],
        how='inner',
        suffixes=('_rMATS', '_LeafCutter')
    )
    return merged

def merge_overlap_genomic(lc_df, rmats_introns_df, max_boundary_drift=100):
    """
    Joins spatial overlaps and filters results using a maximum boundary drift threshold (in bp).
    Enforces strand-matching and prevents short inclusion introns from cross-matching 
    with long exclusion introns.
    """
    # Create PyRanges objects
    pr_lc = pr.PyRanges(lc_df)
    pr_rmats = pr.PyRanges(rmats_introns_df)
    
    # Perform strand-aware spatial join
    overlap_df = pr_rmats.join(pr_lc, suffix='_LeafCutter').df
    
    if overlap_df.empty:
        return overlap_df
        
    # Calculate coordinate differences
    overlap_df['Start_Drift'] = (overlap_df['Start'] - overlap_df['Start_LeafCutter']).abs()
    overlap_df['End_Drift'] = (overlap_df['End'] - overlap_df['End_LeafCutter']).abs()
    
    # Filter: retain matches where BOTH boundaries fall within max_boundary_drift window
    filtered_matches = overlap_df[
        (overlap_df['Start_Drift'] <= max_boundary_drift) & 
        (overlap_df['End_Drift'] <= max_boundary_drift)
    ].copy()
    
    return filtered_matches

# Main Execution Flow
if __name__ == "__main__":
    # leafcutter_file = "data/leafcutter/leafcutter_ds_effect_sizes.csv"
    # rmats_se_file = "data/rmats/SE.csv"
    
    # # Load and standardize
    # lc_data = parse_leafcutter(leafcutter_file)
    # rmats_introns = parse_rmats_se(rmats_se_file)
    
    # # Strategy 1: Exact Intron Boundary Matches
    # exact_matches = merge_exact_introns(lc_data, rmats_introns)
    # exact_matches.to_csv("rmats_leafcutter_exact_matches.csv", index=False)
    
    # # Strategy 2: Drift-Tolerant Overlap Matches (allows up to 30 bp splice site drift)
    # overlap_matches = merge_overlap_genomic(lc_data, rmats_introns, max_boundary_drift=30)
    # overlap_matches.to_csv("rmats_leafcutter_overlap_matches.csv", index=False)

    # Load data
    lc_data = parse_leafcutter("data/leafcutter/leafcutter_ds_effect_sizes.csv")
    rmats_introns = parse_rmats_se("data/rmats/SE.csv")

    print(f"Total rMATS Introns Parsed: {len(rmats_introns)}")
    print(f"Total LeafCutter Introns Parsed: {len(lc_data)}")

    # Step 1: Check Strand compatibility
    print("rMATS Strands:", rmats_introns['Strand'].unique())
    print("LeafCutter Strands:", lc_data['Strand'].unique())

    # Step 2: Unstranded Spatial Overlaps
    pr_lc = pr.PyRanges(lc_data)
    pr_rmats = pr.PyRanges(rmats_introns)

    unstranded_overlaps = pr_rmats.join(pr_lc, suffix='_LeafCutter').df
    print(f"Overlaps found (Unstranded): {len(unstranded_overlaps)}")

    # Step 3: Stranded Spatial Overlaps
    stranded_overlaps = pr_rmats.join(pr_lc, suffix='_LeafCutter').df
    print(f"Overlaps found (Stranded): {len(stranded_overlaps)}")

    # Step 4: After Boundary Drift Filter (30 bp)
    if not stranded_overlaps.empty:
        stranded_overlaps['Start_Drift'] = (stranded_overlaps['Start'] - stranded_overlaps['Start_LeafCutter']).abs()
        stranded_overlaps['End_Drift'] = (stranded_overlaps['End'] - stranded_overlaps['End_LeafCutter']).abs()
        
        filtered = stranded_overlaps[
            (stranded_overlaps['Start_Drift'] <= 30) & 
            (stranded_overlaps['End_Drift'] <= 30)
        ]
        print(f"Matches remaining after 30 bp drift filter: {len(filtered)}")
