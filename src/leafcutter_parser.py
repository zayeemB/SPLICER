

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


