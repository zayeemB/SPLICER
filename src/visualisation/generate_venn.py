"""Venn diagram of junction overlap between LeafCutter, rMATS and MAJIQ.

Run from anywhere:  python src/visualisation/generate_venn.py  (after run.py)

Each tool contributes its filtered, unique junctions (same filters as run.py). Junctions from
different tools are linked when they are on the same strand and within MAX_BOUNDARY_DRIFT bp at both
ends; linked junctions form one "locus", and each locus is placed in the Venn region given by the set
of tools that detected it. The final merged consensus from run.py corresponds to the centre region.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib.pyplot as plt
import pandas as pd
import pyranges as pr
from matplotlib_venn import venn3

from src.config import DATASETS_DIR, MAX_BOUNDARY_DRIFT, PLOTS_DIR
from src.pipeline import load_filtered_tools

MERGED_FILE = DATASETS_DIR / "multitool_summary.csv"
OUTPUT_PATH = PLOTS_DIR / "multitool_venn_diagram.png"
TOOLS = ['LeafCutter', 'rMATS', 'MAJIQ']


def unique_junctions(df):
    """One row per distinct junction (rMATS repeats a junction across events)."""
    cols = ['Chromosome', 'Strand', 'Start', 'End']
    return df[cols].drop_duplicates().reset_index(drop=True)


def link_junctions(tool_junctions):
    """Groups junctions across tools into loci. Returns {frozenset of tools: number of loci}."""
    # Union-find over (tool, row) nodes
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for tool, df in tool_junctions.items():
        for i in df.index:
            find((tool, i))

    names = list(tool_junctions)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            left = tool_junctions[a].rename_axis('idx_a').reset_index()
            right = tool_junctions[b].rename_axis('idx_b').reset_index()
            joined = pr.PyRanges(left).join(pr.PyRanges(right), suffix='_b').df
            if joined.empty:
                continue
            close = joined[
                ((joined['Start'] - joined['Start_b']).abs() <= MAX_BOUNDARY_DRIFT)
                & ((joined['End'] - joined['End_b']).abs() <= MAX_BOUNDARY_DRIFT)
            ]
            for ia, ib in zip(close['idx_a'], close['idx_b']):
                parent[find((a, ia))] = find((b, ib))

    members = {}
    for node in list(parent):
        members.setdefault(find(node), set()).add(node[0])

    counts = {}
    for tools in members.values():
        counts[frozenset(tools)] = counts.get(frozenset(tools), 0) + 1
    return counts


def main():
    lc, rmats, majiq, _ = load_filtered_tools()
    tool_junctions = {
        'LeafCutter': unique_junctions(lc),
        'rMATS': unique_junctions(rmats),
        'MAJIQ': unique_junctions(majiq),
    }
    for tool, df in tool_junctions.items():
        print(f"{tool}: {len(df)} unique filtered junctions")

    counts = link_junctions(tool_junctions)

    def n(*tools):
        return counts.get(frozenset(tools), 0)

    # venn3 order: (A only, B only, A&B, C only, A&C, B&C, A&B&C) with A=LeafCutter, B=rMATS, C=MAJIQ
    subsets = (
        n('LeafCutter'), n('rMATS'), n('LeafCutter', 'rMATS'), n('MAJIQ'),
        n('LeafCutter', 'MAJIQ'), n('rMATS', 'MAJIQ'), n('LeafCutter', 'rMATS', 'MAJIQ'),
    )
    print(f"Venn subsets (loci): {subsets}")

    merged_rows = len(pd.read_csv(MERGED_FILE)) if MERGED_FILE.exists() else None
    shared_by_all = subsets[-1]

    fig, ax = plt.subplots(figsize=(8, 8))
    venn = venn3(subsets=subsets, set_labels=TOOLS, ax=ax)
    for subset_id in ['100', '010', '001', '110', '101', '011', '111']:
        patch = venn.get_patch_by_id(subset_id)
        if patch:
            patch.set_alpha(0.5)

    ax.set_title("Junctions shared across tools (after per-tool filtering)",
                 fontsize=14, fontweight="bold", pad=20)
    caption = f"Shared by all three tools: {shared_by_all} loci"
    if merged_rows is not None:
        caption += f"  |  final merged table: {merged_rows} rows"
    caption += f"\nJunctions linked within {MAX_BOUNDARY_DRIFT} bp at both ends, same strand"
    fig.text(0.5, 0.06, caption, ha='center', fontsize=10)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Venn diagram to '{OUTPUT_PATH}'")


if __name__ == "__main__":
    main()
