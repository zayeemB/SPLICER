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
    
    overlap_df = pr_rmats.join(pr_lc, stranded=True, suffix='_LeafCutter').df
    
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
    leafcutter_file = "leafcutter_effect_sizes.csv"
    
    # Load LeafCutter data
    lc_data = parse_leafcutter(leafcutter_file)
    
    # Load and parse all rMATS event files
    rmats_dfs = []
    
    event_files = {
        'SE': ('SE.MATS.JC.csv', parse_rmats_se),
        'A3SS': ('A3SS.MATS.JC.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A3SS')),
        'A5SS': ('A5SS.MATS.JC.csv', lambda p: parse_rmats_a3ss_a5ss(p, 'A5SS')),
        'RI': ('RI.MATS.JC.csv', parse_rmats_ri),
        'MXE': ('MXE.MATS.JC.csv', parse_rmats_mxe)
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

        # Define your significance thresholds (adjust these numbers as needed)
        pvalue_threshold = 0.05
        fdr_threshold = 0.05
        
        # 1. Filter and save by P-Value
        matches_pvalue = overlap_matches[overlap_matches['PValue'] <= pvalue_threshold].copy()
        matches_pvalue.to_csv("rmats_leafcutter_matches_pvalue_filtered.csv", index=False)
        print(f"Saved {len(matches_pvalue)} matches using P-Value <= {pvalue_threshold}")
        
        # 2. Filter and save by FDR
        matches_fdr = overlap_matches[overlap_matches['FDR'] <= fdr_threshold].copy()
        matches_fdr.to_csv("rmats_leafcutter_matches_fdr_filtered.csv", index=False)
        print(f"Saved {len(matches_fdr)} matches using FDR <= {fdr_threshold}")

        print(f"Total validated matches across all event types: {len(overlap_matches)}")