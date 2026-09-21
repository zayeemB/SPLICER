import pandas as pd
import pyranges as pr

def parse_leafcutter(leafcutter_path):
    """
    Parses LeafCutter CSV output.
    Extracts Chromosome, Start (0-based), End (converted to 0-based), and Strand.
    """
    lc_df = pd.read_csv(leafcutter_path)
    
    split_coords = lc_df['intron'].str.split(':', expand=True)
    lc_df['Chromosome'] = split_coords[0]
    lc_df['Start'] = split_coords[1].astype(int)
    lc_df['End'] = split_coords[2].astype(int) - 1
    lc_df['Strand'] = split_coords[3].str.split('_').str[-1]
    
    if not lc_df['Chromosome'].str.startswith('chr').any():
        lc_df['Chromosome'] = 'chr' + lc_df['Chromosome'].astype(str)
        
    lc_df['Strand'] = lc_df['Strand'].replace({'1': '+', '-1': '-'})
    return lc_df

def standardize_chrom_strand(df):
    """Standardizes chromosome prefix and strand format for rMATS tables."""
    df['Chromosome'] = df['chr']
    if not df['Chromosome'].str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    df['Strand'] = df['strand']
    return df

def parse_rmats_se(path):
    """Parses Skipped Exon (SE) events."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    intron_rows = []
    for _, row in df.iterrows():
        event_id = row['ID']
        gene = row.get('geneSymbol', '')
        chrom, strand = row['Chromosome'], row['Strand']
        diff, fdr, pval = row['IncLevelDifference'], row['FDR'], row['PValue']
        
        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'SE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['upstreamEE'], 
            'End': row['exonStart_0base'], 
            'Intron_Type': 'SE_Upstream', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'SE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['exonEnd'], 
            'End': row['downstreamES'], 
            'Intron_Type': 'SE_Downstream', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'SE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['upstreamEE'], 
            'End': row['downstreamES'], 
            'Intron_Type': 'SE_Exclusion', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

    return pd.DataFrame(intron_rows)

def parse_rmats_a3ss_a5ss(path, event_type):
    """Parses Alternative 3' (A3SS) and Alternative 5' (A5SS) Splice Site events."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    intron_rows = []
    for _, row in df.iterrows():
        event_id = row['ID']
        gene = row.get('geneSymbol', '')
        chrom, strand = row['Chromosome'], row['Strand']
        diff, fdr, pval = row['IncLevelDifference'], row['FDR'], row['PValue']
        
        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': event_type, 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['flankingES'], 
            'End': row['shortES'], 
            'Intron_Type': f'{event_type}_Flanking_Short', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': event_type, 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['flankingES'], 
            'End': row['longExonStart_0base'], 
            'Intron_Type': f'{event_type}_Flanking_Long', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': event_type, 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['shortEE'], 
            'End': row['flankingEE'], 
            'Intron_Type': f'{event_type}_Short_Flanking', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': event_type, 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['longExonEnd'], 
            'End': row['flankingEE'], 
            'Intron_Type': f'{event_type}_Long_Flanking', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

    return pd.DataFrame(intron_rows)

def parse_rmats_ri(path):
    """Parses Retained Intron (RI) events."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    intron_rows = []
    for _, row in df.iterrows():
        event_id = row['ID']
        gene = row.get('geneSymbol', '')
        chrom, strand = row['Chromosome'], row['Strand']
        diff, fdr, pval = row['IncLevelDifference'], row['FDR'], row['PValue']
        
        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'RI', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['upstreamEE'], 
            'End': row['downstreamES'], 
            'Intron_Type': 'RI_SplicedOut', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

    return pd.DataFrame(intron_rows)

def parse_rmats_mxe(path):
    """Parses Mutually Exclusive Exon (MXE) events."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    intron_rows = []
    for _, row in df.iterrows():
        event_id = row['ID']
        gene = row.get('geneSymbol', '')
        chrom, strand = row['Chromosome'], row['Strand']
        diff, fdr, pval = row['IncLevelDifference'], row['FDR'], row['PValue']
        
        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'MXE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['upstreamEE'], 
            'End': row['1stExonStart_0base'], 
            'Intron_Type': 'MXE_Upstream_Exon1', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'MXE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['upstreamEE'], 
            'End': row['2ndExonStart_0base'], 
            'Intron_Type': 'MXE_Upstream_Exon2', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'MXE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['1stExonEnd'], 
            'End': row['downstreamES'], 
            'Intron_Type': 'MXE_Exon1_Downstream', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

        intron_rows.append({
            'rMATS_ID': event_id, 
            'EventType': 'MXE', 
            'geneSymbol': gene, 
            'Chromosome': chrom, 
            'Strand': strand, 
            'Start': row['2ndExonEnd'], 
            'End': row['downstreamES'], 
            'Intron_Type': 'MXE_Exon2_Downstream', 
            'IncLevelDifference': diff, 
            'FDR': fdr,
            'PValue': pval
        })

    return pd.DataFrame(intron_rows)

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

if __name__ == "__main__":
    DATA_DIR = "./data"
    OUTPUT_DIR = "./output"

    leafcutter_file = f"{DATA_DIR}/leafcutter/leafcutter_ds_effect_sizes.csv"
    
    # Load LeafCutter data
    lc_data = parse_leafcutter(leafcutter_file)
    
    # Load and parse all rMATS event files
    rmats_dfs = []
    
    event_files = {
        'SE': (f'{DATA_DIR}/rmats/SE.csv', parse_rmats_se),
        'A3SS': (f'{DATA_DIR}/rmats/A3SS.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A3SS')),
        'A5SS': (f'{DATA_DIR}/rmats/A5SS.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A5SS')),
        'RI': (f'{DATA_DIR}/rmats/RI.csv', parse_rmats_ri),
        'MXE': (f'{DATA_DIR}/rmats/MXE.csv', parse_rmats_mxe)
    }
    
    for event_type, (filename, parser_func) in event_files.items():
        try:
            df_parsed = parser_func(filename)
            rmats_dfs.append(df_parsed)
            print(f"Successfully parsed {len(df_parsed)} introns from {event_type}")
        except FileNotFoundError:
            print(f"Skipping {event_type}: File {filename} not found.")
            
    if rmats_dfs:
        all_rmats_introns = pd.concat(rmats_dfs, ignore_index=True)
        print(f"Total rMATS Introns Across All Event Types: {len(all_rmats_introns)}")
        
        # Run drift-tolerant spatial join against LeafCutter
        overlap_matches = merge_overlap_genomic(lc_data, all_rmats_introns, max_boundary_drift=30)

        # Define your significance and effect size thresholds
        PVALUE_THRESHOLD = 0.01      # Adjust as needed (e.g., 0.01 or 0.05)
        FDR_THRESHOLD = 0.05         # Standard multiple testing threshold
        MIN_DELTA_PSI = 0.10         # Minimum absolute change in inclusion level (10%)

        # 1. P-Value Filtered Dataset
        # Filters for events meeting the p-value cutoff AND showing a real magnitude change
        pvalue_filtered = overlap_matches[
            (overlap_matches['PValue'] < PVALUE_THRESHOLD) & 
            (overlap_matches['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
        ].copy()

        pvalue_filtered.to_csv(f"{OUTPUT_DIR}/rmats_leafcutter_pvalue_filtered.csv", index=False)
        print(f"Saved {len(pvalue_filtered)} rows to 'rmats_leafcutter_pvalue_filtered.csv'")


        # 2. FDR Filtered Dataset
        # Filters for events meeting the False Discovery Rate cutoff AND showing a real magnitude change
        fdr_filtered = overlap_matches[
            (overlap_matches['FDR'] < FDR_THRESHOLD) & 
            (overlap_matches['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
        ].copy()

        fdr_filtered.to_csv(f"{OUTPUT_DIR}/rmats_leafcutter_fdr_filtered.csv", index=False)
        print(f"Saved {len(fdr_filtered)} rows to 'rmats_leafcutter_fdr_filtered.csv'")

        print(f"Total validated matches across all event types: {len(overlap_matches)}")