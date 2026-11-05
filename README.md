# σD Stress Adaptation Models

Decision-tree analysis and growth modeling for studying how variants of the housekeeping sigma factor σD (RpoD) affect *Escherichia coli* adaptation to acid, osmotic, and isobutanol stress.

Associated study: *Natural evolution drives general stress fitness by attenuating the assembly of bacterial house-keeping transcriptional machinery*.

## Models

| Folder | Purpose | Main scripts |
| --- | --- | --- |
| [Decision-tree model](decision_tree_model/) | Relate σD variant features to stress-adaptation phenotypes. | `train.py` (Supplementary Note S1), `train_v2.py`, `train_v3.py` |
| [Growth model](growth_model/) | Compare mutant and WT growth using time-weighted growth-rate differences across normal and stressed conditions. | `calculate_growth_model.py`, `plot_growth_heatmap.py` |

The decision-tree workflow uses single-tree classifiers. The study's GENIE3 random-forest analysis of transcriptional regulatory networks is a separate analysis and is not included here.

## Getting started

Use Python 3.12 and install the requirements for the workflow you need:

```bash
python -m pip install -r growth_model/requirements.txt
python -m pip install -r decision_tree_model/requirements.txt
```

Input formats and example commands are provided in each model folder. Experimental input files must be supplied separately.

## Code and data availability

The analysis code is available in this repository under the [MIT License](LICENSE). Experimental data, model result tables, and result figures are not included.
