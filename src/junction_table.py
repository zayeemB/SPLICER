import numpy as np
import pandas as pd

JUNCTION_KEY = ['Chromosome', 'Strand', 'Start', 'End']


def _join_unique(values):
    return ';'.join(sorted({str(v) for v in values if pd.notna(v) and str(v) != ''}))


def _has_both_signs(series):
    signs = np.sign(series)
    return bool((signs > 0).any() and (signs < 0).any())


def build_junction_table(merged):
    """Collapses the event-level merged table to one row per junction without picking a record.

    A junction is identified by its LeafCutter coordinates. Where a junction belongs to several rMATS
    events or MAJIQ LSVs, their DeltaPSI values are summarised as a min-max range, along with the best
    significance and flags for sign conflicts, instead of choosing one record.
    """
    rows = []
    for key, grp in merged.groupby(JUNCTION_KEY, sort=False):
        lc = grp['IncLevelDifference'].iloc[0]  # one value per LeafCutter junction
        rm, mj = grp['IncLevelDifference_rMATS'], grp['IncLevelDifference_MAJIQ']

        all_signs = np.sign(pd.concat([pd.Series([lc]), rm, mj]))
        rows.append({
            'Chromosome': key[0],
            'Strand': key[1],
            'Start': key[2],
            'End': key[3],
            'Gene': _join_unique(grp['geneSymbol']) or _join_unique(grp['geneSymbol_MAJIQ']),
            'rMATS_EventTypes': _join_unique(grp['EventType']),
            'N_rMATS_Events': grp[['EventType', 'rMATS_ID']].drop_duplicates().shape[0],
            'N_MAJIQ_LSVs': grp['Event_ID'].nunique(),
            'LeafCutter_DeltaPSI': lc,
            'rMATS_DeltaPSI_Min': rm.min(),
            'rMATS_DeltaPSI_Max': rm.max(),
            'MAJIQ_DeltaPSI_Min': mj.min(),
            'MAJIQ_DeltaPSI_Max': mj.max(),
            'rMATS_Best_FDR': grp['FDR'].min(),
            'MAJIQ_Best_ProbabilityChanging': grp['ProbabilityChanging'].max(),
            'rMATS_Sign_Conflict': _has_both_signs(rm),
            'MAJIQ_Sign_Conflict': _has_both_signs(mj),
            'All_Signs_Agree': all_signs.nunique() == 1,
        })
    return pd.DataFrame(rows)
