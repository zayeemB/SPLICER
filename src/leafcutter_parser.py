import pandas as pd

def parse_leafcutter(leafcutter_path):
    """
    Parses LeafCutter effect sizes CSV output.
    Extracts Chromosome, Start (0-based), End (0-based half-open), Strand,
    and maps LeafCutter's DeltaPSI to 'IncLevelDifference' for pipeline compatibility.
    """
    lc_df = pd.read_csv(leafcutter_path)
    
    # Split the standard LeafCutter intron ID column (chrom:start:end:cluster_id)
    split_coords = lc_df['intron'].str.split(':', expand=True)
    lc_df['Chromosome'] = split_coords[0]
    lc_df['Start'] = split_coords[1].astype(int)
    lc_df['End'] = split_coords[2].astype(int) - 1  # Convert 1-based inclusive to 0-based half-open
    lc_df['Strand'] = split_coords[3].str.split('_').str[-1]
    
    # Ensure chromosome format matches (e.g., 'chr1' instead of '1')
    if not lc_df['Chromosome'].str.startswith('chr').any():
        lc_df['Chromosome'] = 'chr' + lc_df['Chromosome'].astype(str)
        
    lc_df['Strand'] = lc_df['Strand'].replace({'1': '+', '-1': '-'})
    
    # Dynamically find and map LeafCutter's DeltaPSI column to 'IncLevelDifference'
    delta_psi_col = next((c for c in lc_df.columns if c.lower() in ['deltapsi', 'delta_psi']), None)
    if delta_psi_col:
        lc_df['IncLevelDifference'] = lc_df[delta_psi_col]
    else:
        raise ValueError(f"LeafCutter file is missing a DeltaPSI column. Available columns: {list(lc_df.columns)}")
        
    lc_df['Tool_Source'] = 'LeafCutter'
    return lc_df

