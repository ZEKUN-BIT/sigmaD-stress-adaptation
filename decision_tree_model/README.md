# Decision-tree model

Analyze relationships between σD variant features and strain adaptation using `DecisionTreeClassifier`. The main entry point, `train.py`, follows the nine-feature, threshold-1 workflow embedded in the original manuscript's Supplementary Note S1 and supplied separately in Note S1 and `BIO/train.py`.

The manuscript's three trees refer to acid, osmotic, and isobutanol stress: run the main workflow separately with labels for each condition. `train_v2.py` and `train_v3.py` retain additional feature-subset variants from the BIO archive; they do not represent those three stress conditions.

## Input

Provide two Excel workbooks: training data with an ID in the first column, features in the middle, and a numeric source label in the last column; and candidate data containing only features in the same order. The first row contains headers.

| Script | Input features, in order | Default label threshold |
| --- | --- | --- |
| `train.py` | Original residue, variant residue, upstream sequence, downstream sequence, domain, position, core-polymerase affinity change, DNA affinity change, stability change | 1 |
| `train_v2.py` | Original residue, variant residue, domain, upstream sequence, downstream sequence, core-polymerase affinity change, DNA affinity change, stability change; the two sequence columns are omitted during fitting | 2 |
| `train_v3.py` | Original residue, variant residue, domain, position, core-polymerase affinity change, DNA affinity change, stability change | 2 |

Source labels strictly above the threshold become class 1; all others become class 0. The supplied ordinal encoding and 80:20 split (`random_state=0`) are retained. The classifier retains the source default parameters, including an unset tree random seed.

The threshold applies to numeric source labels before binary conversion, rather than to the 0/1 classes described in the manuscript's Methods. Already binarized 0/1 labels would all become class 0 at the Note S1 threshold of 1.

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

The scripts cover the classifier and these exports. Feature-importance plots, partial-dependence plots, and the separate GENIE3 network analysis described in the manuscript are not implemented here.

The original manuscript reports Python 3.10.13, scikit-learn 1.3.2, and pandas 2.1.3. This packaged version was checked with Python 3.12.3 and the versions pinned in `requirements.txt`; exact agreement with historical results is not guaranteed across software versions or runs with the unset tree random seed.
