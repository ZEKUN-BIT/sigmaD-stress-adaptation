#!/usr/bin/env python3
"""Calculate growth advantages and thresholds without plotting.

Run: python calculate_growth_model.py --input source_data.csv
Dependencies: numpy, pandas
All rates and rate differences are in min^-1.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REFERENCE = "p\u03c3-MG"
STRAINS = [f"p\u03c3-{i}" for i in range(1, 43)]
CONDITIONS = ["Normal", "Acid", "NaCl", "IBOH"]


def load_rates(path):
    data = pd.read_csv(path, float_precision="round_trip")
    required = {"strain", "condition", "mu_input", "input_unit"}
    if not required.issubset(data.columns):
        raise ValueError(f"Required columns: {sorted(required)}")
    if set(data["input_unit"]) != {"min^-1"}:
        raise ValueError("Input rates must use min^-1.")

    # Retain strain IDs and condition codes, dropping descriptive suffixes.
    data["strain"] = data["strain"].str.extract(
        r"^([^-]+-(?:MG|\d+))(?=$|[^A-Za-z0-9])", expand=False
    )
    data["condition"] = data["condition"].str.extract(
        r"^(Normal|Acid|NaCl|IBOH)(?=$|[^A-Za-z0-9])", expand=False
    )
    if data[["strain", "condition"]].isna().any().any():
        raise ValueError("Unrecognized strain or condition label.")
    if data.duplicated(["strain", "condition"]).any():
        raise ValueError("Expected one value per strain and condition.")
    if (set(data["strain"]) != {REFERENCE, *STRAINS}
            or set(data["condition"]) != set(CONDITIONS)):
        raise ValueError("Expected WT, mutants 1-42, and four conditions.")

    rates = data.pivot(index="strain", columns="condition", values="mu_input")
    rates = rates.loc[[REFERENCE, *STRAINS], CONDITIONS]
    values = rates.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("All 172 rates must be finite and nonnegative.")
    return rates.astype(float)


def calculate(rates):
    delta = rates.loc[STRAINS] - rates.loc[REFERENCE]
    p_values = np.linspace(0, 1, 101)
    rows, scenarios = [], []
    for stress in CONDITIONS[1:]:
        for strain in STRAINS:
            a = delta.loc[strain, "Normal"]
            b = delta.loc[strain, stress]
            p_star = np.nan
            if a * b < 0:
                p_star = -a / (b - a)
                if b > 0:
                    kind = "normal_cost_stress_benefit"
                    advantage = "p > p_star"
                else:
                    kind = "normal_benefit_stress_cost"
                    advantage = "p < p_star"
            elif a > 0 and b > 0:
                kind, advantage = "beneficial_both", "all p"
            elif a < 0 and b < 0:
                kind, advantage = "deleterious_both", "none"
            else:
                raise ValueError(
                    "Endpoint ties or neutral rows require additional "
                    f"classification rules: {strain}, {stress}."
                )
            rows.append({
                "strain": strain,
                "stress": stress,
                "delta_normal_min_inv": a,
                "delta_stress_min_inv": b,
                "s_25_min_inv": 0.75 * a + 0.25 * b,
                "s_50_min_inv": 0.5 * (a + b),
                "s_75_min_inv": 0.25 * a + 0.75 * b,
                "classification": kind,
                "p_star": p_star,
                "p_star_percent": 100 * p_star,
                "advantage": advantage,
            })
            scenarios.extend({
                "strain": strain, "stress": stress,
                "p_stress": p, "s_min_inv": (1 - p) * a + p * b,
            } for p in p_values)
    return pd.DataFrame(rows), pd.DataFrame(scenarios)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Path to a user-supplied input CSV.")
    parser.add_argument("--output", type=Path, default=Path("results/calculation_results"))
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input CSV not found: {args.input}")
    results, scenarios = calculate(load_rates(args.input))
    counts = results.groupby(["stress", "classification"]).size().unstack(fill_value=0)
    args.output.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output / "model_results.csv", index=False, encoding="utf-8-sig")
    scenarios.to_csv(args.output / "scenario_values.csv", index=False, encoding="utf-8-sig")
    counts.to_csv(args.output / "classification_counts.csv", encoding="utf-8-sig")
    print(f"Saved {len(results)} comparisons and {len(scenarios)} scenarios.")
    print(f"Output directory: {args.output.resolve()}")


if __name__ == "__main__":
    main()
