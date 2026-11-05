# Environment-weighted growth model

Compare mutant and WT growth as the fraction of time under stress changes. Each stress condition is analyzed separately. The calculation uses time-weighted mutant-minus-WT growth-rate differences and returns point estimates and transition thresholds.

## Model quantities

For the assumptions, symbol definitions, and step-by-step derivation from the population-growth equation, see [Model definition and derivation](MODEL.md).

For each mutant and one selected stress condition, define the growth-rate differences:

$$
a = \mu_{\mathrm{mut,normal}} - \mu_{\mathrm{WT,normal}}, \qquad
b = \mu_{\mathrm{mut,stress}} - \mu_{\mathrm{WT,stress}}.
$$

Let $p$ be the fraction of total time spent under that stress, with the remaining fraction $1-p$ spent under normal conditions. The model calculates:

$$
s(p) = (1-p)a + pb, \qquad 0 \leq p \leq 1.
$$

The quantities $a$, $b$, and $s(p)$ all have units of min⁻¹. Positive $s(p)$ indicates a mutant growth advantage, negative $s(p)$ indicates a WT advantage, and zero indicates equal time-weighted growth rates. This is a growth-rate difference, not a percentage change or the logarithm of a measured growth rate. Here, $p$ is an environmental time fraction, not a statistical p-value.

The term *relative log-growth rate* follows from exponential growth. If $R(t)=N_{\mathrm{mut}}(t)/N_{\mathrm{WT}}(t)$ is the population ratio, then:

$$
\ln\left(\frac{R(T)}{R(0)}\right) = T\,s(p).
$$

Here, $T$ is total time in minutes and $\ln$ is the natural logarithm. At fixed $p$, $s(p)$ is independent of $T$, while the accumulated population ratio depends on $T$. The scripts calculate rates and thresholds; they do not require $T$ or output population sizes.

With equal starting populations, $R(0)=1$ and $R(T)=\exp[T s(p)]$. With unequal starting populations, $R(T)=R(0)\exp[T s(p)]$: positive $s(p)$ means that the ratio increases, without necessarily becoming greater than 1.

This interpretation treats each supplied maximum growth rate as a constant rate within its condition. It does not model lag phases, resource limitation, or changes in growth rate with population density. The results are point estimates; confidence intervals and statistical significance are not calculated.

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

## Output files and fields

| File | Contents |
| --- | --- |
| `model_results.csv` | 126 rows: one summary for each of 42 mutants under each of 3 stresses. |
| `scenario_values.csv` | 12,726 rows: $s(p)$ for every mutant–stress pair at 101 fractions, from 0 to 1 in steps of 0.01. |
| `classification_counts.csv` | 3 rows: mutant counts by classification for each stress. Each row sums to 42. Columns represent classes observed in the dataset; an absent class within a stress has count 0. |
| `threshold_summary.csv` | Plotting output: the `strain`, `stress`, `delta_normal_min_inv`, `delta_stress_min_inv`, `classification`, `p_star`, and `advantage` columns from the summary calculation. |
| `heatmap_values.csv` | Plotting output: the same columns and values as `scenario_values.csv` when supplied the same input. |

The fields in `model_results.csv` are:

| Field | Meaning |
| --- | --- |
| `strain` | Mutant identifier; WT is the reference for every comparison. |
| `stress` | Selected stress: `Acid`, `NaCl`, or `IBOH`. Stresses are evaluated separately. |
| `delta_normal_min_inv` | $a$: mutant minus WT growth rate under normal conditions, in min⁻¹. |
| `delta_stress_min_inv` | $b$: mutant minus WT growth rate under the selected stress, in min⁻¹. |
| `s_25_min_inv`, `s_50_min_inv`, `s_75_min_inv` | $s(p)$ at stress time fractions of 0.25, 0.50, and 0.75, respectively, in min⁻¹. |
| `classification` | One of the four sign-based categories below. |
| `p_star` | Stress time fraction where $s(p)=0$, when an interior crossover exists. Blank when both rate differences have the same sign. |
| `p_star_percent` | `100 * p_star`, expressed as a percentage; blank when `p_star` is blank. |
| `advantage` | The range of $p$ where the mutant has a strictly positive advantage: `p > p_star`, `p < p_star`, `all p`, or `none`. |

`scenario_values.csv` and `heatmap_values.csv` contain `strain`, `stress`, `p_stress`, and `s_min_inv`. The last two fields are the dimensionless fraction $p$ and its corresponding $s(p)$ in min⁻¹. For example, `p_stress = 0.25` means 25% of time under stress and 75% under normal conditions.

## Thresholds and classifications

When $a$ and $b$ have opposite signs, the crossover is calculated analytically:

$$
p^* = \frac{-a}{b-a}.
$$

This threshold lies strictly between 0 and 1 and is not rounded to the 0.01 scenario grid. At $p=p^*$, neither strain has a modeled growth advantage.

| `classification` | Rate differences | Mutant advantage |
| --- | --- | --- |
| `normal_cost_stress_benefit` | $a<0$, $b>0$ | `p > p_star` |
| `normal_benefit_stress_cost` | $a>0$, $b<0$ | `p < p_star` |
| `beneficial_both` | $a>0$, $b>0$ | `all p`, including 0 and 1 |
| `deleterious_both` | $a<0$, $b<0$ | `none` |

A blank threshold means there is no crossover within the modeled range; it does not indicate missing input data. Cases with $a=0$ or $b=0$ are rejected by the current classification rules.

For an illustrative example using artificial values, let $a=-0.001$ and $b=0.002$ min⁻¹. Then $p^*=1/3$, so the mutant has an advantage when stress occupies more than approximately 33.3% of the time. At $p=0.5$, $s(p)=0.0005$ min⁻¹. With equal starting populations and $T=1000$ minutes, the exponential-growth interpretation gives $R(T)=\exp(0.5)\approx1.65$.

## Reading the heatmap

The plotting script saves `growth_advantage_heatmap_Arial_original_style` in PNG, PDF, and SVG formats. Rows are mutants, columns are stress time fractions shown as percentages, and the three panels represent separate stress conditions.

Red indicates positive $s(p)$ and blue indicates negative $s(p)$, using a shared scale centered on zero. The colorbar displays values in units of 10⁻³ min⁻¹: a displayed value of 1 corresponds to `s_min_inv = 0.001`. CSV values remain in min⁻¹ and are not multiplied by 1000. No normalization is applied to individual rows.

The `>` and `<` markers locate $p^*$ and indicate whether the mutant advantage occurs above or below that threshold. In the `All p` strip, `+` means mutant advantage throughout the range and `−` means WT advantage throughout. Colors and markers describe model point estimates, not statistical significance.

## Tests

```bash
cd growth_model
python -m unittest discover -s tests -v
```

Tests use artificial inputs only.
