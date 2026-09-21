import pandas as pd
import pyranges as pr

def standardize_chrom_strand(df):
    """Standardizes chromosome prefix and strand format for rMATS tables."""
    df['Chromosome'] = df['chr']
    if not df['Chromosome'].str.startswith('chr').any():
        df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    df['Strand'] = df['strand']
    return df

def merge_overlap_genomic(df1, df2, max_boundary_drift=30, suffix='_right'):
    """
    Performs strand-aware PyRanges spatial join with a maximum boundary drift threshold,
    using a dynamic suffix to prevent column collisions during sequential joins.
    """
    pr_1 = pr.PyRanges(df1)
    pr_2 = pr.PyRanges(df2)
    
    # Perform spatial join using the dynamic suffix
    overlap_df = pr_1.join(pr_2, suffix=suffix).df
    
    if overlap_df.empty:
        return overlap_df
        
    # Dynamically compute boundary drift using the specified suffix
    overlap_df['Start_Drift'] = (overlap_df['Start'] - overlap_df[f'Start{suffix}']).abs()
    overlap_df['End_Drift'] = (overlap_df['End'] - overlap_df[f'End{suffix}']).abs()
    
    filtered_matches = overlap_df[
        (overlap_df['Start_Drift'] <= max_boundary_drift) & 
        (overlap_df['End_Drift'] <= max_boundary_drift)
    ].copy()
    
    return filtered_matches

def load_gene_chromosome_map(gtf_path):
    """
    Loads a reference GTF file (e.g., GENCODE basic CHR GTF) and builds 
    a dictionary mapping Ensembl gene IDs to their reference chromosomes.
    """
    print(f"Loading reference GTF from {gtf_path}...")
    gtf = pr.read_gtf(gtf_path)
    df = gtf.df
    
    # Filter strictly for 'gene' feature lines to keep things clean and efficient
    df_genes = df[df['Feature'] == 'gene'] if 'Feature' in df.columns else df
    
    # Build the dictionary mapping gene_id -> Chromosome
    if 'gene_id' in df_genes.columns and 'Chromosome' in df_genes.columns:
        # Strip version numbers from Ensembl IDs if present (e.g., ENSG00000123456.5 -> ENSG00000123456)
        clean_ids = df_genes['gene_id'].astype(str).str.split('.').str[0]
        gene_map = dict(zip(clean_ids, df_genes['Chromosome']))
    else:
        raise ValueError("GTF file is missing expected 'gene_id' or 'Chromosome' columns.")
        
    print(f"Successfully built chromosome map for {len(gene_map)} genes.")
    return gene_map