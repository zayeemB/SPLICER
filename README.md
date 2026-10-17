# SPLICER

**Splice-junction & Pathway Locus Integration for Coordinate-matched Events in RNA-Seq**

SPLICER finds differential splicing events that several RNA-seq tools agree on. It reads the output
of **LeafCutter**, **rMATS** (SE, A3SS, A5SS, RI, MXE) and **MAJIQ**, puts them on a common footing,
and reports the junctions that all three tools call as changing between hypoxia and normoxia.

## How it works

1. **Parse.** Each tool's output becomes one table of junctions (chromosome, strand, 0-based
   half-open start/end, DeltaPSI, significance). rMATS events are split into their individual
   junctions, and each junction's DeltaPSI sign is set by the isoform it supports.
2. **Same contrast.** DeltaPSI is **Hypoxia - Normoxia** for every tool. rMATS (group 1 = hypoxia)
   already is; LeafCutter and MAJIQ report Norm - Hyp and are flipped.
3. **Filter each tool before merging.**
   - rMATS: FDR < 0.05 and |DeltaPSI| >= 0.10
   - MAJIQ: probability changing >= 0.90 and |DeltaPSI| >= 0.10
   - LeafCutter: |DeltaPSI| >= 0.10
4. **Match junctions.** Junctions are matched by coordinates (same strand, within 2 bp at both ends,
   closest match only): LeafCutter with rMATS first, then with MAJIQ.
5. **Report.** Sign agreement between tools, junctions that passed the filters but were not merged,
   and the read support of the unmerged rMATS junctions (raw and depth-normalised).

All thresholds and paths are in `src/config.py`.

## Setup

Python 3 with `pandas`, `numpy`, `scipy`, `pyranges`, `matplotlib` and `matplotlib-venn`. The project
includes a virtual environment in `.env/`:

```bash
source .env/bin/activate
```

## Input data

Put the tool outputs under `data/`:

| Path | Content |
|---|---|
| `data/leafcutter/leafcutter_ds_effect_sizes.csv` | LeafCutter effect sizes (needs `intron` and `deltapsi` columns) |
| `data/majiq/tsv_f.csv` | MAJIQ Voila TSV/CSV (needs `seqid`, `junctions_coords`, `mean_dpsi_per_lsv_junction`, `probability_changing`) |
| `data/rmats/{SE,A3SS,A5SS,RI,MXE}.csv` | rMATS event files, with group 1 = hypoxia |

A missing file is skipped with a message. To use different file names, edit `src/config.py`.

## Usage

Both scripts work from any directory. Run the pipeline first, then the plots:

```bash
python run.py
python make_plots.py
```

`make_plots.py` stops with a message if `run.py` hasn't produced its tables yet.

## Output

**`output/datasets/`** (written by `run.py`)

| File | Content |
|---|---|
| `multitool_streamlined_filtered.csv` | Full merged table, one row per matched event, all columns from each tool |
| `multitool_summary.csv` | Condensed view: coordinates, chromosome, gene, DeltaPSI and significance per tool |
| `multitool_junctions.csv` | One row per junction, with DeltaPSI ranges and sign-conflict flags |
| `multitool_report.txt` | Filter counts, sign agreement, unmerged junctions, rMATS read support |
| `sign_disagreements.csv` | Merged rows where the three tools disagree on the DeltaPSI sign |

**`output/plots/`** (written by `make_plots.py`)

| File | Content |
|---|---|
| `multitool_venn_diagram.png` | Junction overlap between the three tools after filtering |
| `top_genes_deltapsi_comparison.png` | Top 30 genes by splicing change, DeltaPSI per tool |
| `event_type_comparison.png` | Merge rate, sign agreement and effect size by rMATS event type |

## Project layout

```
run.py                  pipeline: load, filter, merge, save tables and report
make_plots.py           runs every plot script
src/
  config.py             thresholds and paths
  pipeline.py           loading and filtering of the three tools
  parsers/              LeafCutter, rMATS and MAJIQ parsers
  utils.py              coordinate matching
  analysis.py           unmerged-junction and read-support analysis
  junction_table.py     per-junction summary table
  report.py             text report
  visualisation/        one script per plot (each has a main())
```

To add a plot, give a new module in `src/visualisation/` a `main()` function and list it in
`make_plots.py`.

## Notes

- rMATS is the only tool whose output here has per-junction read counts, so the read-support
  analysis covers rMATS only. Library sizes are approximated from the total event reads.
- A junction can belong to several rMATS events or MAJIQ LSVs, so the merged table can have several
  rows per junction. `multitool_junctions.csv` has one row per junction.
