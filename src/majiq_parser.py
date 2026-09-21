import pandas as pd

def parse_majiq_voila(voila_path, gene_to_chrom_map=None):
    """
    Parses MAJIQ/Voila output using the updated schema.
    - Automatically detects file delimiter (CSV/TSV).
    - Extracts chromosomes natively from the 'seqid' column.
    - Maps 'mean_dpsi_per_lsv_junction' and 'probability_changing'.
    """
    # Automatically sniff delimiter (handles comma-separated CSV or tab-separated TSV)
    df = pd.read_csv(voila_path, sep=None, engine='python')
    
    intron_rows = []
    skipped_count = 0
    
    for _, row in df.iterrows():
        raw_gene = str(row.get('gene_id', ''))
        clean_gene = raw_gene.replace('gene:', '').split('.')[0].strip()
        
        # 1. Resolve chromosome using the explicit 'seqid' column from the new schema
        raw_chrom = row.get('seqid', None)
        if pd.notna(raw_chrom):
            raw_chrom_str = str(raw_chrom)
            chrom = 'chr' + raw_chrom_str if not raw_chrom_str.startswith('chr') else raw_chrom_str
        elif gene_to_chrom_map is not None:
            # Fallback to GTF map if seqid is missing
            chrom = gene_to_chrom_map.get(clean_gene, None)
        else:
            chrom = None
            
        if not chrom:
            skipped_count += 1
            continue  # Skip if chromosome cannot be resolved
            
        lsv_id = str(row.get('lsv_id', ''))
        lsv_type = str(row.get('lsv_type', ''))
        strand = str(row.get('strand', '+'))
        
        # 2. Extract semicolon-delimited arrays using the updated column headers
        junc_str = str(row.get('junctions_coords', ''))
        dpsi_str = str(row.get('mean_dpsi_per_lsv_junction', ''))
        prob_str = str(row.get('probability_changing', ''))
        
        if not junc_str or pd.isna(junc_str):
            continue
            
        junc_coords_list = junc_str.split(';')
        dpsi_vals = dpsi_str.split(';') if not pd.isna(dpsi_str) else []
        prob_vals = prob_str.split(';') if not pd.isna(prob_str) else []
        
        for idx, junc_pair in enumerate(junc_coords_list):
            if '-' not in junc_pair:
                continue
            
            parts = junc_pair.split('-')
            start_1base = int(parts[0])
            end_1base = int(parts[1])
            
            # Convert 1-based inclusive coordinates to 0-based half-open format
            start_0base = start_1base - 1
            end_0base = end_1base
            
            dpsi = float(dpsi_vals[idx]) if idx < len(dpsi_vals) else 0.0
            prob = float(prob_vals[idx]) if idx < len(prob_vals) else 0.0
            
            intron_rows.append({
                'Tool_Source': 'MAJIQ',
                'Event_ID': lsv_id,
                'geneSymbol': clean_gene,
                'Chromosome': chrom,
                'Strand': strand,
                'Start': start_0base,
                'End': end_0base,
                'Intron_Type': f'MAJIQ_{lsv_type}',
                'IncLevelDifference': dpsi,          # Mapped from mean_dpsi_per_lsv_junction
                'ProbabilityChanging': prob          # Mapped from probability_changing
            })
            
    if skipped_count > 0:
        print(f"Note: Skipped {skipped_count} MAJIQ rows due to missing chromosome/seqid information.")
        
    return pd.DataFrame(intron_rows)