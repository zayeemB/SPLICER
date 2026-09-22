import pandas as pd

def parse_leafcutter(leafcutter_path, flip_sign=False):
    """
    Parses LeafCutter effect sizes CSV output and standardizes to 
    0-based half-open coordinates [Start, End).
    - Includes a `flip_sign` flag to invert DeltaPSI if contrast groups are reversed.
    """
    lc_df = pd.read_csv(leafcutter_path)
    
    # Split the standard LeafCutter intron ID column (chrom:start:end:cluster_id)
    split_coords = lc_df['intron'].str.split(':', expand=True)
    lc_df['Chromosome'] = split_coords[0]
    
    # CORRECTED 0-based half-open conversion:
    lc_df['Start'] = split_coords[1].astype(int) - 1  # Subtract 1 from 1-based start
    lc_df['End'] = split_coords[2].astype(int)        # Keep 1-based end as exclusive boundary
    
    lc_df['Strand'] = split_coords[3].str.split('_').str[-1]
    
    # Ensure chromosome format matches ('chr' prefix)
    if not lc_df['Chromosome'].str.startswith('chr').any():
        lc_df['Chromosome'] = 'chr' + lc_df['Chromosome'].astype(str)
        
    lc_df['Strand'] = lc_df['Strand'].replace({'1': '+', '-1': '-'})
    
    # Map DeltaPSI to IncLevelDifference
    delta_psi_col = next((c for c in lc_df.columns if c.lower() in ['deltapsi', 'delta_psi']), None)
    if delta_psi_col:
        dpsi = lc_df[delta_psi_col]
        
        # Apply sign flip if contrast groups are inverted
        if flip_sign:
            dpsi = -dpsi
            
        lc_df['IncLevelDifference'] = dpsi
    else:
        raise ValueError(f"LeafCutter file is missing a DeltaPSI column. Available columns: {list(lc_df.columns)}")
        
    lc_df['Tool_Source'] = 'LeafCutter'
    return lc_df