"""Compares the rMATS event types (SE, A3SS, A5SS, RI, MXE) across the merged data.

Run from anywhere:  python src/visualisation/event_type_comparison.py  (after run.py)

Panels:
  A. Junctions passing the rMATS filters vs. how many reached the merged data (and the % that did)
  B. Percentage of merged rows where all three tools agree on the DeltaPSI sign
  C. |DeltaPSI| per tool, by event type
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from src.config import DATASETS_DIR, PLOTS_DIR
from src.pipeline import load_filtered_tools
from src.report import TOOLS, sign_agreement

MERGED_FILE = DATASETS_DIR / "multitool_summary.csv"
OUTPUT_PATH = PLOTS_DIR / "event_type_comparison.png"
EVENT_ORDER = ['SE', 'A3SS', 'A5SS', 'RI', 'MXE']
COLORS = {'LeafCutter': '#e8605d', 'rMATS': '#5aa469', 'MAJIQ': '#5b6ee8'}
GREY, DARK = '#bdbdbd', '#4a4a4a'


def main():
    merged = pd.read_csv(MERGED_FILE)
    _, rmats, _, _ = load_filtered_tools()

    events = [e for e in EVENT_ORDER if e in set(merged['EventType']) | set(rmats['EventType'])]
    x = np.arange(len(events))

    # Distinct junctions per event type, before and after merging
    key = ['EventType', 'Chromosome', 'Strand']
    filtered_n = rmats.drop_duplicates(key + ['Start', 'End']).groupby('EventType').size()
    merged_n = merged.drop_duplicates(key + ['rMATS_Start', 'rMATS_End']).groupby('EventType').size()
    filtered_n = filtered_n.reindex(events, fill_value=0)
    merged_n = merged_n.reindex(events, fill_value=0)

    merged['All_Agree'] = sign_agreement(merged)
    agree_pct = merged.groupby('EventType')['All_Agree'].mean().reindex(events) * 100
    rows_n = merged.groupby('EventType').size().reindex(events, fill_value=0)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # A. filtered vs merged junction counts
    ax = axes[0]
    w = 0.38
    ax.bar(x - w / 2, filtered_n, w, color=GREY, label='Pass rMATS filters')
    ax.bar(x + w / 2, merged_n, w, color=COLORS['rMATS'], label='In merged data')
    for i, (f, m) in enumerate(zip(filtered_n, merged_n)):
        ax.text(i + w / 2, m, f"{100 * m / f:.0f}%" if f else "", ha='center', va='bottom', fontsize=9)
    ax.set_yscale('log')
    ax.set_ylim(bottom=10)
    ax.set_xticks(x, events)
    ax.set_ylabel("Unique rMATS junctions (log scale)")
    ax.set_title("A. Junctions reaching the merged data\n(label = % of filtered)", fontweight='bold')
    ax.legend()

    # B. sign agreement
    ax = axes[1]
    ax.bar(x, agree_pct, color=DARK)
    for i, (p, n) in enumerate(zip(agree_pct, rows_n)):
        ax.text(i, p + 0.5, f"{p:.1f}%\n(n={n})", ha='center', va='bottom', fontsize=9)
    ax.set_ylim(80, 106)
    ax.set_xticks(x, events)
    ax.set_ylabel("Merged rows where all 3 tools agree on sign (%)\n(axis starts at 80%)")
    ax.set_title("B. DeltaPSI sign agreement", fontweight='bold')

    # C. |DeltaPSI| per tool by event type
    ax = axes[2]
    width = 0.26
    for j, tool in enumerate(TOOLS):
        data = [merged.loc[merged['EventType'] == e, f'{tool}_DeltaPSI'].abs() for e in events]
        bp = ax.boxplot(data, positions=x + (j - 1) * width, widths=width * 0.9, patch_artist=True,
                        showfliers=False, medianprops={'color': 'black'})
        for box in bp['boxes']:
            box.set_facecolor(COLORS[tool])
            box.set_alpha(0.8)
    ax.set_xticks(x, events)
    ax.set_ylabel("|DeltaPSI| (Hyp - Norm)")
    ax.set_title("C. Effect size by tool", fontweight='bold')
    ax.legend(handles=[Patch(facecolor=COLORS[t], alpha=0.8, label=t) for t in TOOLS])

    fig.suptitle("Comparison by rMATS event type", fontsize=15, fontweight='bold')
    fig.tight_layout()

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved event-type comparison to '{OUTPUT_PATH}'")


if __name__ == "__main__":
    main()
