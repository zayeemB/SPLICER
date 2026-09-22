import pandas as pd
import pyranges as pr


def determine_tool_sign_flip(parser_func, tool_path, rmats_path, **kwargs):
    """Automatically determines if a tool's DeltaPSI signs need to be inverted

    by checking the Spearman correlation of overlapping events against rMATS.

    Parameters:
    - parser_func (callable): The parsing function (e.g., parse_leafcutter or parse_majiq_voila)
    - tool_path (str): File path to the tool's output file.
    - rmats_path (str): File path to rMATS SE output.
    - **kwargs: Any extra arguments the parser might need.

    Returns:
    - bool: True if signs need to be flipped, False otherwise.
    """
    from src.parsers.rmats_parser import parse_rmats_se

    tool_name = parser_func.__name__
    print(
        f"Running automated sign-alignment check between {tool_name} and rMATS..."
    )

    # 1. Load temporary dataset with sign-flipping explicitly disabled
    temp_tool = parser_func(tool_path, flip_sign=False, **kwargs)
    temp_rmats = parse_rmats_se(rmats_path)

    if temp_tool.empty or temp_rmats.empty:
        print(
            f"-> Warning: One of the datasets is empty for {tool_name}. Defaulting"
            " to False."
        )
        return False

    # 2. Keep the metadata columns inside PyRanges from the start
    pr_tool = pr.PyRanges(temp_tool)
    pr_rmats = pr.PyRanges(temp_rmats)

    # PyRanges join automatically carries over all metadata columns 
    # (rMATS columns will get a '_b' suffix)
    overlaps_pr = pr_tool.join(pr_rmats)
    merged = overlaps_pr.df.drop_duplicates(subset=["Chromosome", "Start", "End"])

    if merged.empty:
        print(f"-> Warning: No coordinate overlaps found for {tool_name}. Defaulting to False.")
        return False

    # 3. Directly rename the columns for correlation check
    # 'IncLevelDifference' belongs to temp_tool, 'IncLevelDifference_b' belongs to temp_rmats
    merged = merged.rename(columns={
        "IncLevelDifference": "Tool_dPSI",
        "IncLevelDifference_b": "rMATS_dPSI"
    })

    if len(merged) < 5:
        print(
            f"-> Warning: Only {len(merged)} overlaps found. Not enough data for"
            f" stable correlation for {tool_name}. Defaulting to False."
        )
        return False

    # 4. Calculate Spearman correlation across overlapping DeltaPSIs
    corr = merged["Tool_dPSI"].corr(merged["rMATS_dPSI"], method="spearman")
    print(f"-> Empirical Spearman Correlation ({tool_name} vs rMATS): {corr:.3f}")

    # 5. Return True if inversion is required, False otherwise
    if corr < 0:
        print(
            "-> Negative correlation detected! Contrasts are inverted. (Returning"
            " True)"
        )
        return True
    else:
        print(
            "-> Positive correlation detected! Contrasts are aligned. (Returning"
            " False)"
        )
        return False


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