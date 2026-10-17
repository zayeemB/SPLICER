"""Shared configuration and loading/filtering of the three tools' outputs.

Used by run.py and the visualisation scripts so they all apply identical filters.
All tools report DeltaPSI as Hypoxia - Normoxia after parsing: rMATS (group 1 = hypoxia) is natively
Hyp - Norm; LeafCutter and MAJIQ are natively Norm - Hyp and are flipped.
"""
from pathlib import Path

import pandas as pd

from .parsers.leafcutter_parser import parse_leafcutter
from .parsers.majiq_parser import parse_majiq_voila
from .parsers.rmats_parser import (
    parse_rmats_se,
    parse_rmats_a3ss_a5ss,
    parse_rmats_mxe,
    parse_rmats_ri,
)

# Statistical and effect size thresholds
RMATS_FDR_THRESHOLD = 0.05       # rMATS q-value (FDR) cutoff
MAJIQ_PROB_THRESHOLD = 0.90      # MAJIQ Probability of Changing cutoff (>= 90% confidence)
MIN_DELTA_PSI = 0.10             # Minimum absolute magnitude change (10%)
MAX_BOUNDARY_DRIFT = 2           # Max bp difference at each junction end (tools differ by ~1 bp by convention)

# All paths are resolved relative to the project root, so scripts work from any working directory.
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "output"

LEAFCUTTER_FILE = DATA_DIR / "leafcutter" / "leafcutter_ds_effect_sizes.csv"
MAJIQ_FILE = DATA_DIR / "majiq" / "tsv_f.csv"
RMATS_DIR = DATA_DIR / "rmats"

RMATS_EVENT_FILES = {
    'SE': (RMATS_DIR / "SE.csv", parse_rmats_se),
    'A3SS': (RMATS_DIR / "A3SS.csv", lambda p: parse_rmats_a3ss_a5ss(p, 'A3SS')),
    'A5SS': (RMATS_DIR / "A5SS.csv", lambda p: parse_rmats_a3ss_a5ss(p, 'A5SS')),
    'RI': (RMATS_DIR / "RI.csv", parse_rmats_ri),
    'MXE': (RMATS_DIR / "MXE.csv", parse_rmats_mxe),
}


def load_filtered_tools():
    """Parses each tool and applies its significance / effect-size filter.

    Returns (lc_data, rmats_data, majiq_data, filter_counts), where filter_counts maps a tool name
    to (rows before filtering, rows after filtering, description of the filter). A tool whose input
    file is missing comes back as an empty DataFrame.
    """
    filter_counts = {}

    # LeafCutter
    try:
        lc_data = parse_leafcutter(LEAFCUTTER_FILE, flip_sign=True)
        print(f"Successfully parsed {len(lc_data)} introns from LeafCutter.")
        n_before = len(lc_data)
        lc_data = lc_data[lc_data['IncLevelDifference'].abs() >= MIN_DELTA_PSI]
        filter_counts['LeafCutter'] = (n_before, len(lc_data), f"abs(DeltaPSI) >= {MIN_DELTA_PSI}")
        print(f"  -> {len(lc_data)} pass |DeltaPSI| >= {MIN_DELTA_PSI}")
    except FileNotFoundError:
        lc_data = pd.DataFrame()
        print(f"Skipping LeafCutter: File {LEAFCUTTER_FILE} not found.")

    # MAJIQ
    try:
        majiq_data = parse_majiq_voila(MAJIQ_FILE, flip_sign=True)
        print(f"Successfully parsed and mapped {len(majiq_data)} flat introns from MAJIQ.")
        n_before = len(majiq_data)
        majiq_data = majiq_data[
            (majiq_data['ProbabilityChanging'] >= MAJIQ_PROB_THRESHOLD)
            & (majiq_data['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
        ]
        filter_counts['MAJIQ'] = (
            n_before, len(majiq_data),
            f"probability >= {MAJIQ_PROB_THRESHOLD}, abs(DeltaPSI) >= {MIN_DELTA_PSI}")
        print(f"  -> {len(majiq_data)} pass probability >= {MAJIQ_PROB_THRESHOLD} and |DeltaPSI| >= {MIN_DELTA_PSI}")
    except FileNotFoundError:
        majiq_data = pd.DataFrame()
        print(f"Skipping MAJIQ: File {MAJIQ_FILE} not found.")

    # rMATS (all event types)
    rmats_dfs = []
    for event_type, (filename, parser_func) in RMATS_EVENT_FILES.items():
        try:
            df_parsed = parser_func(filename)
            rmats_dfs.append(df_parsed)
            print(f"Successfully parsed {len(df_parsed)} introns from rMATS {event_type}")
        except FileNotFoundError:
            print(f"Skipping rMATS {event_type}: File {filename} not found.")

    if rmats_dfs:
        rmats_data = pd.concat(rmats_dfs, ignore_index=True)
        print(f"Total rMATS Introns Across All Event Types: {len(rmats_data)}")
        n_before = len(rmats_data)
        rmats_data = rmats_data[
            (rmats_data['FDR'] < RMATS_FDR_THRESHOLD)
            & (rmats_data['IncLevelDifference'].abs() >= MIN_DELTA_PSI)
        ]
        filter_counts['rMATS'] = (
            n_before, len(rmats_data),
            f"FDR < {RMATS_FDR_THRESHOLD}, abs(DeltaPSI) >= {MIN_DELTA_PSI}")
        print(f"  -> {len(rmats_data)} pass FDR < {RMATS_FDR_THRESHOLD} and |DeltaPSI| >= {MIN_DELTA_PSI}")
    else:
        rmats_data = pd.DataFrame()

    return lc_data, rmats_data, majiq_data, filter_counts
