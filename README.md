# Code for σD variant analysis

This repository contains the two analysis workflows used in the study of σD variants and strain adaptation to stress.

## Models

| Folder | Purpose | Main scripts |
| --- | --- | --- |
| [Growth model](growth_model/) | Compare mutant and WT growth across different stress time fractions. | `calculate_growth_model.py`, `plot_growth_heatmap.py` |
| [Decision-tree model](decision_tree_model/) | Analyze variant features associated with stress adaptation, following Supplementary Note S1. | `train.py`, `train_v2.py`, `train_v3.py` |

```text
.
├── growth_model/          # Growth calculations and heatmap code
├── decision_tree_model/   # Decision-tree training and prediction code
├── LICENSE
└── README.md
```

## Getting started

Use Python 3.12 and install the requirements for the workflow you need:

```bash
python -m pip install -r growth_model/requirements.txt
python -m pip install -r decision_tree_model/requirements.txt
```

Input formats and example commands are provided in each model folder. Experimental input files must be supplied separately.

## Code and data availability

The analysis code is available in this repository under the [MIT License](LICENSE). Experimental data, model result tables, and result figures are not included.
