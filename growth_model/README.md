# Environment-weighted growth model

Compare mutant and WT growth as the fraction of time under stress changes. Each stress condition is analyzed separately. The calculation uses time-weighted mutant-minus-WT growth-rate differences and returns point estimates and transition thresholds.

This is the revised growth model developed during manuscript revision. It replaces the original Note S2 population-size simulation and is not a reproduction of the original Fig. 5j.

## Input

A UTF-8 CSV with `strain`, `condition`, `mu_input`, and `input_unit` columns. Supply one value for each of WT `pσ-MG` and mutants `pσ-1` to `pσ-42` under `Normal`, `Acid`, `NaCl`, and `IBOH` (172 rows). Rates must be finite and nonnegative, with `input_unit` set to `min^-1`.

Missing or duplicate comparisons and endpoint ties with WT are rejected. Uncertainty is not estimated.

## Run

From the repository root:

```bash
python -m pip install -r growth_model/requirements.txt
python growth_model/calculate_growth_model.py --input /path/to/rates.csv --output outputs/growth
python growth_model/plot_growth_heatmap.py --input /path/to/rates.csv --output outputs/growth_figures
```

The calculation saves comparison, scenario, and classification-count CSVs. Plotting saves the heatmap as PNG, PDF, and SVG, together with its underlying values. Arial must be installed separately for plotting.

## Tests

```bash
cd growth_model
python -m unittest discover -s tests -v
```

Tests use artificial inputs only.
