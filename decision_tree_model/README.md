# Decision-tree model

Analyze relationships between σD variant features and strain adaptation using the single-tree classifier supplied in Supplementary Note S1. The two additional entry points retain the feature subsets supplied in the BIO code archive. All three use `DecisionTreeClassifier`.

## Input

Provide two Excel workbooks: training data with an ID in the first column, features in the middle, and a numeric source label in the last column; and candidate data containing only features in the same order. The first row contains headers.

| Script | Input features, in order | Default label threshold |
| --- | --- | --- |
| `train.py` | Original residue, variant residue, upstream sequence, downstream sequence, domain, position, core-polymerase affinity change, DNA affinity change, stability change | 1 |
| `train_v2.py` | Original residue, variant residue, domain, upstream sequence, downstream sequence, core-polymerase affinity change, DNA affinity change, stability change; the two sequence columns are omitted during fitting | 2 |
| `train_v3.py` | Original residue, variant residue, domain, position, core-polymerase affinity change, DNA affinity change, stability change | 2 |

Source labels strictly above the threshold become class 1; all others become class 0. The supplied ordinal encoding and 80:20 split (`random_state=0`) are retained. The classifier retains the source default parameters, including an unset tree random seed.

## Run

From the repository root:

```bash
python -m pip install -r decision_tree_model/requirements.txt
python decision_tree_model/train.py --train /path/to/training.xlsx --predict /path/to/candidates.xlsx --output outputs/decision_tree
```

Use `train_v2.py` or `train_v3.py` with the matching input layout for the corresponding feature subset. Use `--threshold` to override its default.

Outputs are held-out predictions, candidate predictions, per-class precision/recall/F1 CSVs, and `decision_tree.dot`. The CSV export includes both the source label and binary class. If Graphviz is installed, render the exported tree with:

```bash
dot -Tpdf outputs/decision_tree/decision_tree.dot -o outputs/decision_tree/decision_tree.pdf
```

The implementation consolidates the supplied scripts, makes input/output paths explicit, and corrects the validation CSV header count. No input workbooks or generated outputs are distributed.
