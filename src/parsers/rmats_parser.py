import pandas as pd

from ..utils import standardize_chrom_strand

# Contrast convention for every parser in this package: IncLevelDifference = Hypoxia - Normoxia.
# rMATS was run with group 1 = hypoxia, so its IncLevelDifference (IncLevel1 - IncLevel2) is already
# Hyp - Norm *for the inclusion isoform*. Junctions that belong to the alternative (skipping /
# short / spliced-out) isoform move the opposite way, so their sign is inverted here.

OUTPUT_COLS = [
    'rMATS_ID', 'EventType', 'geneSymbol', 'Chromosome', 'Strand', 'Start', 'End',
    'Intron_Type', 'IncLevelDifference', 'FDR', 'PValue',
]


def _load(path, coord_cols):
    """Reads an rMATS table, standardises chrom/strand and drops rows with unusable coordinates."""
    df = pd.read_csv(path)
    df = standardize_chrom_strand(df)
    for col in coord_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=coord_cols + ['IncLevelDifference']).copy()
    df[coord_cols] = df[coord_cols].astype(int)
    return df


def _junctions(df, event_type, specs):
    """Builds one output row per (event, junction).

    specs: list of (Intron_Type, start, end, sign) where start/end are column names (or Series) in
    df and sign is +1 for inclusion-isoform junctions and -1 for alternative-isoform junctions.
    Introns are 0-based half-open: Start = end of the upstream exon, End = 0-based start of the next.
    """
    frames = []
    for intron_type, start, end, sign in specs:
        out = pd.DataFrame({
            'rMATS_ID': df['ID'],
            'EventType': event_type,
            'geneSymbol': df.get('geneSymbol', ''),
            'Chromosome': df['Chromosome'],
            'Strand': df['Strand'],
            'Start': df[start] if isinstance(start, str) else start,
            'End': df[end] if isinstance(end, str) else end,
            'Intron_Type': intron_type,
            'IncLevelDifference': sign * df['IncLevelDifference'],
            'FDR': df['FDR'],
            'PValue': df['PValue'],
        })
        frames.append(out)
    return pd.concat(frames, ignore_index=True)[OUTPUT_COLS]


def parse_rmats_se(path):
    """Parses Skipped Exon (SE) events. Inclusion = exon included."""
    df = _load(path, ['upstreamEE', 'exonStart_0base', 'exonEnd', 'downstreamES'])
    return _junctions(df, 'SE', [
        ('SE_Upstream', 'upstreamEE', 'exonStart_0base', +1),
        ('SE_Downstream', 'exonEnd', 'downstreamES', +1),
        ('SE_Exclusion', 'upstreamEE', 'downstreamES', -1),
    ])


def parse_rmats_a3ss_a5ss(path, event_type):
    """Parses A3SS / A5SS events. Inclusion = long exon, alternative = short exon.

    The flanking exon lies upstream (genomically left) of the alternative exon for A3SS on '+' and
    A5SS on '-', and to the right otherwise, so the junction coordinates depend on that side.
    """
    cols = ['flankingES', 'flankingEE', 'shortES', 'shortEE', 'longExonStart_0base', 'longExonEnd']
    df = _load(path, cols)

    left = df['flankingEE'] <= df['longExonStart_0base']  # flanking exon is genomically left
    flank_start = df['flankingEE'].where(left, df['shortEE'])
    short_end = df['shortES'].where(left, df['flankingES'])
    long_start = df['flankingEE'].where(left, df['longExonEnd'])
    long_end = df['longExonStart_0base'].where(left, df['flankingES'])

    return _junctions(df, event_type, [
        (f'{event_type}_Long', long_start, long_end, +1),
        (f'{event_type}_Short', flank_start, short_end, -1),
    ])


def parse_rmats_ri(path):
    """Parses Retained Intron (RI) events. Inclusion = intron retained, so the spliced junction is -1."""
    df = _load(path, ['upstreamEE', 'downstreamES'])
    return _junctions(df, 'RI', [
        ('RI_SplicedOut', 'upstreamEE', 'downstreamES', -1),
    ])


def parse_rmats_mxe(path):
    """Parses Mutually Exclusive Exon (MXE) events.

    Per the rMATS docs the inclusion isoform is the 1st exon on '+' and the 2nd exon on '-'.
    """
    cols = ['upstreamEE', '1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',
            'downstreamES']
    df = _load(path, cols)

    sign_1st = df['Strand'].map({'+': 1, '-': -1}).fillna(1)
    sign_2nd = -sign_1st

    frames = []
    for intron_type, start, end, sign in [
        ('MXE_Upstream_Exon1', 'upstreamEE', '1stExonStart_0base', sign_1st),
        ('MXE_Exon1_Downstream', '1stExonEnd', 'downstreamES', sign_1st),
        ('MXE_Upstream_Exon2', 'upstreamEE', '2ndExonStart_0base', sign_2nd),
        ('MXE_Exon2_Downstream', '2ndExonEnd', 'downstreamES', sign_2nd),
    ]:
        part = _junctions(df, 'MXE', [(intron_type, start, end, 1)])
        part['IncLevelDifference'] = part['IncLevelDifference'] * sign.values
        frames.append(part)
    return pd.concat(frames, ignore_index=True)[OUTPUT_COLS]
