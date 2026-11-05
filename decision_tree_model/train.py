#!/usr/bin/env python3
"""Train the decision-tree workflow supplied in Supplementary Note S1.

Run: python train.py --train training.xlsx --predict candidates.xlsx
The v2 and v3 entry points retain the feature subsets in the BIO archive.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, precision_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder
from sklearn.tree import DecisionTreeClassifier, export_graphviz


FEATURE_SETS = {
    "full": {"input_features": 9, "categorical": 5, "numeric_start": 5, "threshold": 1},
    "v2": {"input_features": 8, "categorical": 3, "numeric_start": 5, "threshold": 2},
    "v3": {"input_features": 7, "categorical": 3, "numeric_start": 3, "threshold": 2},
}


def run_analysis(train_path, prediction_path, output, feature_set="full", threshold=None):
    """Fit the source classifier and export held-out and candidate predictions."""
    settings = FEATURE_SETS[feature_set]
    if threshold is None:
        threshold = settings["threshold"]
    training = pd.read_excel(train_path, header=0)
    candidates = pd.read_excel(prediction_path, header=0)
    n_features = settings["input_features"]
    if training.shape[1] != n_features + 2 or candidates.shape[1] != n_features:
        raise ValueError(
            f"The {feature_set} workflow expects {n_features + 2} training columns "
            f"(ID, {n_features} features, label) and {n_features} candidate columns."
        )

    raw_x = training.iloc[:, 1:-1].to_numpy()
    raw_y = training.iloc[:, -1].to_numpy()
    raw_candidates = candidates.to_numpy()
    n_categorical = settings["categorical"]
    numeric_start = settings["numeric_start"]
    selected = list(range(n_categorical)) + list(range(numeric_start, n_features))
    feature_names = training.columns[1:-1][selected].tolist()

    # Preserve the supplied workflow: fit the encoder before the held-out split.
    encoder = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1)
    encoded = encoder.fit_transform(raw_x[:, :n_categorical])
    x = np.concatenate((encoded, raw_x[:, numeric_start:]), axis=1)
    encoded_candidates = encoder.transform(raw_candidates[:, :n_categorical])
    candidate_x = np.concatenate((encoded_candidates, raw_candidates[:, numeric_start:]), axis=1)

    x_train, x_test, source_train, source_test = train_test_split(
        x, raw_y, test_size=0.2, random_state=0
    )
    y_train = np.where(source_train > threshold, 1, 0)
    y_test = np.where(source_test > threshold, 1, 0)

    # Keep the original single-tree estimator and its default parameters.
    classifier = DecisionTreeClassifier()
    classifier.fit(x_train, y_train)
    heldout_predictions = classifier.predict(x_test)
    candidate_predictions = classifier.predict(candidate_x)

    decoded_test = np.concatenate(
        (encoder.inverse_transform(x_test[:, :n_categorical]), x_test[:, n_categorical:]),
        axis=1,
    )
    validation = pd.DataFrame(decoded_test, columns=feature_names)
    validation["source_label"] = source_test
    validation["class"] = y_test
    validation["prediction"] = heldout_predictions
    predictions = candidates.iloc[:, selected].copy()
    predictions.columns = feature_names
    predictions["prediction"] = candidate_predictions
    metrics = pd.DataFrame({
        "class": np.unique(np.concatenate((y_test, heldout_predictions))),
        "precision": precision_score(y_test, heldout_predictions, average=None, zero_division=0),
        "recall": recall_score(y_test, heldout_predictions, average=None, zero_division=0),
        "f1": f1_score(y_test, heldout_predictions, average=None, zero_division=0),
    })

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    validation.to_csv(output / "validation_predictions.csv", index=False, encoding="utf-8-sig")
    predictions.to_csv(output / "candidate_predictions.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(output / "classification_metrics.csv", index=False, encoding="utf-8-sig")
    export_graphviz(
        classifier, out_file=str(output / "decision_tree.dot"),
        feature_names=[str(name) for name in feature_names],
    )
    return classifier


def main(feature_set="full"):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train", type=Path, required=True, help="Labeled training workbook.")
    parser.add_argument("--predict", type=Path, required=True, help="Unlabeled candidate workbook.")
    parser.add_argument("--output", type=Path, default=Path("results/decision_tree"))
    parser.add_argument(
        "--threshold", type=float, default=FEATURE_SETS[feature_set]["threshold"],
        help="Assign class 1 only when the source label is strictly above this threshold.",
    )
    args = parser.parse_args()
    for path in [args.train, args.predict]:
        if not path.is_file():
            parser.error(f"Input workbook not found: {path}")
    run_analysis(args.train, args.predict, args.output, feature_set, args.threshold)
    print(f"Saved predictions, classification metrics, and decision_tree.dot to {args.output.resolve()}")


if __name__ == "__main__":
    main()
