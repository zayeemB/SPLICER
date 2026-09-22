import pandas as pd
import pyranges as pr

from ..utils import standardize_chrom_strand

def parse_rmats_se(path):
    """Parses Skipped Exon (SE) events with explicit integer coordinate casting."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    
    # Ensure chromosome format matches LeafCutter/MAJIQ ('chr' prefix)
    if 'Chromosome' in df.columns and not df['Chromosome'].astype(str).str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
        
    # Explicitly cast coordinate columns to integers
    coord_cols = ['upstreamEE', 'exonStart_0base', 'exonEnd', 'downstreamES']
    for col in coord_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

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
    
    if 'Chromosome' in df.columns and not df['Chromosome'].astype(str).str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
        
    coord_cols = ['flankingES', 'shortES', 'longExonStart_0base', 'shortEE', 'longExonEnd', 'flankingEE']
    for col in coord_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

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
    
    if 'Chromosome' in df.columns and not df['Chromosome'].astype(str).str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
        
    coord_cols = ['upstreamEE', 'downstreamES']
    for col in coord_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

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
    
    if 'Chromosome' in df.columns and not df['Chromosome'].astype(str).str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
        
    coord_cols = ['upstreamEE', '1stExonStart_0base', '2ndExonStart_0base', '1stExonEnd', '2ndExonEnd', 'downstreamES']
    for col in coord_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

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