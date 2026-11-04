# Environment-weighted growth model

Python scripts for comparing mutant and wild-type growth when normal and stress conditions alternate. A mutation can improve growth under stress while reducing growth under normal conditions. This model quantifies how the predicted advantage changes with the fraction of time spent under stress.

Each stress condition is paired with the normal condition and analyzed separately. The scripts calculate relative log-growth rates and transition thresholds, and optionally draw a heatmap.

## Code and data availability

This is a **code-only release**. Experimental measurements from the three experiments, input tables, model result tables, and study figures are not included. Users must supply an input CSV to run the analysis on experimental data. No dataset accession or data download link is provided in this repository.

The tests generate artificial values in memory and temporary files. These values are unrelated to the experimental measurements.

## Model definition

For mutant strain $i$ and the wild-type reference (WT), let $\mu_{i,n}$ and $\mu_{i,s}$ denote growth rates under normal and stress conditions. All input rates are in **min⁻¹**. Within each environment, the model assumes exponential growth:

$$
\frac{dN_i}{dt}=\mu_{i,e}N_i,\qquad e\in\{n,s\}.
$$

Let $T_n$ and $T_s$ be the cumulative durations in the two environments, with $T=T_n+T_s>0$. Define the actual stress time fraction as

$$
p=\frac{T_s}{T},\qquad 0\leq p\leq 1.
$$

Integrating growth over the two environments gives

$$
\ln\frac{N_i(T)}{N_i(0)}=\mu_{i,n}T_n+\mu_{i,s}T_s.
$$

Define the mutant-minus-WT rate differences

$$
a_i=\mu_{i,n}-\mu_{\mathrm{WT},n},\qquad
b_i=\mu_{i,s}-\mu_{\mathrm{WT},s}.
$$

The environment-weighted relative log-growth rate is

$$
s_i(p)=(1-p)a_i+pb_i
=\frac{1}{T}\ln\left[\frac{N_i(T)/N_i(0)}{N_{\mathrm{WT}}(T)/N_{\mathrm{WT}}(0)}\right].
$$

A positive $s_i(p)$ predicts mutant advantage; a negative value predicts WT advantage; zero indicates equal growth. The unit of $s_i$ is min⁻¹. The input growth rates are used directly, **without a logarithmic transformation**. The logarithm in the definition concerns population expansion.

For equal initial population sizes,

$$
\frac{N_i(T)}{N_{\mathrm{WT}}(T)}=\exp\bigl[T s_i(p)\bigr].
$$

At a fixed actual time fraction $p$, $s_i(p)$ is independent of total duration, while the accumulated log population ratio $T s_i(p)$ depends on duration. Durations must be expressed in minutes when rates are in min⁻¹. If an alternating schedule ends partway through a cycle, $p$ must reflect the actual cumulative durations, including that partial cycle.

When $a_i$ and $b_i$ have opposite signs, the internal transition threshold is obtained by solving $s_i(p)=0$:

$$
p_i^*=-\frac{a_i}{b_i-a_i}.
$$

| Normal difference | Stress difference | Classification | Mutant advantage |
| --- | --- | --- | --- |
| $a_i>0$ | $b_i>0$ | `beneficial_both` | All $p$ |
| $a_i<0$ | $b_i>0$ | `normal_cost_stress_benefit` | $p>p_i^*$ |
| $a_i>0$ | $b_i<0$ | `normal_benefit_stress_cost` | $p<p_i^*$ |
| $a_i<0$ | $b_i<0$ | `deleterious_both` | No $p$ |

At $p=p_i^*$, the model predicts equal growth. This implementation rejects endpoint ties ($a_i=0$ or $b_i=0$), including fully neutral comparisons, because their classification and plotting symbols are not defined in this version.

### Assumptions and interpretation

- Supplied maximum growth rates are treated as constant rates within each environment. Switching is instantaneous, with no lag, adaptation, or environmental memory.
- The model omits carrying capacity, resource depletion, density effects, death, and interactions between strains. Population ratios compare modeled growth trajectories; they are not a mechanistic simulation of competition in mixed culture.
- Under these assumptions, only cumulative time in each environment affects the prediction; switching order and frequency have no additional effect.
- The scripts return point estimates. They do not fit raw growth curves, estimate confidence intervals, or perform significance tests. If replicate means are supplied, predictions are based on those means. Replicate-level measurements and a suitable analysis design are needed to assess uncertainty.

## Installation

Tested on Linux with Python **3.12.3**, NumPy **2.4.3**, pandas **2.3.3**, and Matplotlib **3.10.8**.

```bash
git clone https://github.com/ZEKUN-BIT/environment-weighted-growth-model.git
cd environment-weighted-growth-model
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Calculation uses NumPy and pandas. Plotting additionally requires Matplotlib and a separately installed **Arial** font. The plotting script checks for Arial and stops if it is unavailable. Font files are not distributed with this code.

## Input format

Supply a UTF-8 CSV with the following columns:

| Column | Meaning |
| --- | --- |
| `strain` | `pσ-MG` for WT; `pσ-1` through `pσ-42` for mutants |
| `condition` | Exactly one of `Normal`, `Acid`, `NaCl`, or `IBOH` |
| `mu_input` | A finite, nonnegative growth rate |
| `input_unit` | The literal string `min^-1` |

The expected layout contains **172 rows**: one value for each of 43 strains under each of four conditions. Missing combinations and duplicate strain-condition pairs are rejected. The scripts do not average replicate rows automatically. Descriptive suffixes after recognized strain or condition identifiers are removed before validation; use the canonical identifiers above for new inputs.

No measured input values are embedded in the source code. Run each experiment with its own input file and a distinct output directory.

## Run the calculation

```bash
python calculate_growth_model.py \
  --input /path/to/your_input.csv \
  --output results/my_analysis
```

The calculation uses an exact analytical expression, so no numerical ODE solver is needed. For each mutant and stress condition, it evaluates $s_i(p)$ on 101 equally spaced points from 0 to 1 and calculates $p_i^*$ analytically when an internal threshold exists.

| Output | Contents |
| --- | --- |
| `model_results.csv` | Normal and stress rate differences; $s_i$ at 25%, 50%, and 75% stress time; classification; threshold; advantage direction |
| `scenario_values.csv` | $s_i(p)$ for every mutant, stress condition, and sampled $p$ |
| `classification_counts.csv` | Number of comparisons in each classification, grouped by stress condition |

Rate columns ending in `_min_inv` are in min⁻¹. `p_star` and `p_stress` are fractions; `p_star_percent` is a percentage. A blank threshold means no internal transition exists for that comparison.

## Draw the heatmap

```bash
python plot_growth_heatmap.py \
  --input /path/to/your_input.csv \
  --output figures/my_analysis
```

The plotting script calls the same calculation functions. It creates one panel per stress condition, with mutants as rows and stress time fraction as columns. All panels share a symmetric `RdBu_r` color scale centered at zero, with displayed values in 10⁻³ min⁻¹. Threshold symbols point toward the region of mutant advantage; the `All p` column marks comparisons with the same advantage throughout the range.

The figure uses Arial, 9.5 pt strain labels and ticks, 12 pt titles, and no panel letters. It is exported at 183 × 210 mm as a 600 dpi PNG and as PDF/SVG with editable text. Output names are `growth_advantage_heatmap_Arial_original_style` with the respective extension. The script also saves `threshold_summary.csv` and `heatmap_values.csv`.

Generated outputs and common data/figure formats are excluded by `.gitignore`.

## Verification

```bash
python -m unittest discover -s tests -v
```

The synthetic tests check advantage directions and thresholds, endpoint values, equivalence to piecewise exponential growth, independence from duration at a fixed time fraction, input validation, and agreement between calculation and plotting values. No experimental data are required.

## Suggested manuscript Code availability statement

> The Python code implementing the environment-weighted growth model and heatmap visualization is available at https://github.com/ZEKUN-BIT/environment-weighted-growth-model (version v0.1.0) under the MIT license. Experimental input data and generated results are not included in the code repository.

This statement describes code access only. Data availability should be reported separately according to the manuscript's actual data access arrangements. Cite the version used for the final analysis; no archival DOI is assigned by this repository.

## License

The code and documentation are released under the [MIT License](LICENSE).
