"""Model verification using artificial values, never experimental data."""

from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from calculate_growth_model import CONDITIONS, REFERENCE, STRAINS, calculate, load_rates
from plot_growth_heatmap import prepare_data


def synthetic_rates():
    """Generate all four sign combinations with a known threshold of 1/3."""
    reference = np.array([0.01, 0.008, 0.006, 0.004])
    differences = [
        [-0.001, 0.002, 0.002, 0.002],
        [0.001, -0.002, -0.002, -0.002],
        [0.001, 0.002, 0.002, 0.002],
        [-0.001, -0.002, -0.002, -0.002],
    ]
    values = [reference]
    values.extend(reference + differences[i % 4] for i in range(len(STRAINS)))
    return pd.DataFrame(values, index=[REFERENCE, *STRAINS], columns=CONDITIONS)


def input_table(rates):
    data = rates.rename_axis("strain").reset_index().melt(
        id_vars="strain", var_name="condition", value_name="mu_input"
    )
    data["input_unit"] = "min^-1"
    return data


class GrowthModelTests(unittest.TestCase):
    def setUp(self):
        self.rates = synthetic_rates()
        self.results, self.scenarios = calculate(self.rates)

    def test_thresholds_and_advantage_directions(self):
        acid = self.results[self.results.stress == "Acid"].set_index("strain")
        expected = [
            ("normal_cost_stress_benefit", "p > p_star"),
            ("normal_benefit_stress_cost", "p < p_star"),
            ("beneficial_both", "all p"),
            ("deleterious_both", "none"),
        ]
        for strain, (kind, advantage) in zip(STRAINS[:4], expected):
            row = acid.loc[strain]
            self.assertEqual(row.classification, kind)
            self.assertEqual(row.advantage, advantage)
        for strain in STRAINS[:2]:
            row = acid.loc[strain]
            self.assertAlmostEqual(row.p_star, 1 / 3)
            self.assertAlmostEqual((1 - row.p_star) * row.delta_normal_min_inv
                                   + row.p_star * row.delta_stress_min_inv, 0)
        self.assertTrue(acid.loc[STRAINS[2:4], "p_star"].isna().all())

    def test_grid_endpoints_and_quartiles(self):
        self.assertEqual(len(self.results), 42 * 3)
        self.assertEqual(len(self.scenarios), 42 * 3 * 101)
        for stress in CONDITIONS[1:]:
            curves = self.scenarios[self.scenarios.stress == stress]
            for fraction, condition in [(0, "Normal"), (1, stress)]:
                observed = curves[curves.p_stress == fraction].set_index("strain")
                expected = self.rates.loc[STRAINS, condition] - self.rates.loc[REFERENCE, condition]
                np.testing.assert_allclose(observed.loc[STRAINS, "s_min_inv"], expected)
        row = self.results.iloc[0]
        np.testing.assert_allclose(
            [row.s_25_min_inv, row.s_50_min_inv, row.s_75_min_inv],
            [-0.00025, 0.0005, 0.00125], atol=1e-17,
        )

    def test_piecewise_population_growth_and_time_scaling(self):
        mutant, wt = self.rates.loc[STRAINS[0]], self.rates.loc[REFERENCE]
        schedule = [("Normal", 15), ("Acid", 30), ("Normal", 30), ("Acid", 45)]
        row = self.results.iloc[0]
        duration = sum(minutes for _, minutes in schedule)
        p = sum(minutes for env, minutes in schedule if env == "Acid") / duration
        s = (1 - p) * row.delta_normal_min_inv + p * row.delta_stress_min_inv
        for multiplier in [1, 2]:
            for order in [schedule, list(reversed(schedule))]:
                n_mutant, n_wt = 1.0, 1.0
                for environment, minutes in order:
                    n_mutant *= np.exp(mutant[environment] * minutes * multiplier)
                    n_wt *= np.exp(wt[environment] * minutes * multiplier)
                log_ratio = np.log(n_mutant / n_wt)
                self.assertAlmostEqual(log_ratio, duration * multiplier * s)
                self.assertAlmostEqual(log_ratio / (duration * multiplier), s)

    def test_input_round_trip_and_descriptive_suffixes(self):
        data = input_table(self.rates)
        data["strain"] += " (description)"
        data["condition"] += " (description)"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.csv"
            data.to_csv(path, index=False)
            loaded = load_rates(path)
        np.testing.assert_array_equal(loaded.to_numpy(), self.rates.to_numpy())
        self.assertEqual(loaded.index.tolist(), [REFERENCE, *STRAINS])

    def test_invalid_inputs_are_rejected(self):
        data = input_table(self.rates)
        cases = {
            "duplicate": pd.concat([data, data.iloc[[0]]], ignore_index=True),
            "missing": data.iloc[1:].copy(),
            "unit": data.assign(input_unit="h^-1"),
            "unknown_condition": data.replace({"condition": {"Acid": "Other"}}),
            "negative": data.assign(mu_input=-0.001),
            "nonfinite": data.assign(mu_input=np.inf),
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.csv"
            for label, invalid in cases.items():
                with self.subTest(case=label):
                    invalid.to_csv(path, index=False)
                    with self.assertRaises(ValueError):
                        load_rates(path)

    def test_endpoint_ties_are_explicitly_rejected(self):
        for condition in ["Normal", "Acid"]:
            with self.subTest(condition=condition):
                rates = self.rates.copy()
                rates.loc[STRAINS[0], condition] = rates.loc[REFERENCE, condition]
                with self.assertRaisesRegex(ValueError, "Endpoint ties"):
                    calculate(rates)

    def test_plot_uses_calculation_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.csv"
            input_table(self.rates).to_csv(path, index=False)
            p, matrices, summary, curves = prepare_data(path)
        pd.testing.assert_frame_equal(curves, self.scenarios)
        pd.testing.assert_frame_equal(summary, self.results[summary.columns])
        for stress, matrix in matrices.items():
            selected = self.scenarios[self.scenarios.stress == stress]
            self.assertEqual(matrix.shape, (42, 101))
            for index, strain in enumerate(STRAINS):
                expected = selected[selected.strain == strain].set_index("p_stress")
                np.testing.assert_array_equal(matrix[index], expected.loc[p, "s_min_inv"])


if __name__ == "__main__":
    unittest.main()
