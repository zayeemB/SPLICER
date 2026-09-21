import pandas as pd
import pyranges as pr

def parse_majiq(voila_path):
    """
    Parses MAJIQ/Voila TSV output, flattening multi-junction LSVs 
    into individual intron rows for pipeline integration.
    """
    df = pd.read_csv(voila_path, sep='\t')
    
    intron_rows = []
    for _, row in df.iterrows():
        gene = str(row.get('gene_id', row.get('Gene name', '')))
        lsv_id = str(row.get('lsv_id', ''))
        lsv_type = str(row.get('lsv_type', ''))
        
        # Extract chromosome and strand 
        chrom = str(row.get('chrom', row.get('Chromosome', row.get('Chr', ''))))
        strand = str(row.get('strand', row.get('Strand', '')))
        
        # Dynamically find columns starting with 'junctions' or containing stats
        junc_col = next((c for c in df.columns if c.startswith('junctions')), None)
        dpsi_col = next((c for c in df.columns if 'mean_dpsi' in c), None)
        prob_col = next((c for c in df.columns if 'probability' in c), None)
        
        if not junc_col or pd.isna(row[junc_col]):
            continue
            
        junc_strs = str(row[junc_col]).split(':')
        dpsi_vals = str(row[dpsi_col]).split(':') if dpsi_col and not pd.isna(row[dpsi_col]) else [0.0] * len(junc_strs)
        prob_vals = str(row[prob_col]).split(':') if prob_col and not pd.isna(row[prob_col]) else [1.0] * len(junc_strs)
        
        for idx, junc_str in enumerate(junc_strs):
            if '-' not in junc_str:
                continue
            parts = junc_str.split('-')
            start_1base = int(parts[0])
            end_1base = int(parts[1])
            
            # Convert 1-based inclusive coordinates to 0-based half-open format
            start_0base = start_1base - 1
            end_0base = end_1base
            
            dpsi = float(dpsi_vals[idx]) if idx < len(dpsi_vals) else 0.0
            prob = float(prob_vals[idx]) if idx < len(prob_vals) else 1.0
            
            formatted_chrom = 'chr' + chrom if not str(chrom).startswith('chr') else chrom
            
            intron_rows.append({
                'Tool_Source': 'MAJIQ',
                'Event_ID': lsv_id,
                'geneSymbol': gene,
                'Chromosome': formatted_chrom,
                'Strand': strand,
                'Start': start_0base,
                'End': end_0base,
                'Intron_Type': f'MAJIQ_{lsv_type}',
                'IncLevelDifference': dpsi,       # Mapped as Delta-PSI
                'ProbabilityChanging': prob,
                'PValue': 1.0 - prob             # Approximate p-value equivalent for filtering compatibility
            })
            
    return pd.DataFrame(intron_rows)