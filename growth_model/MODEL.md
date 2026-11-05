# Model definition and derivation

This model compares one mutant with WT across normal conditions and one selected stress condition. Acid, NaCl, and IBOH are evaluated separately. The reported quantity is the difference between their time-weighted growth rates, in min⁻¹.

## 1. Parameters and assumptions

| Symbol | Definition | Unit or range |
| --- | --- | --- |
| $i$ | Strain: the selected mutant (`mut`) or WT | — |
| $N_i(t)$ | Modeled population size of strain $i$ at time $t$, with $N_i(0)>0$ | Cell number |
| $\mu_{i,n}$ | Growth rate of strain $i$ under normal conditions | min⁻¹ |
| $\mu_{i,s}$ | Growth rate of strain $i$ under the selected stress | min⁻¹ |
| $T_n$, $T_s$ | Total time spent in normal and stressed conditions, respectively | Minutes; nonnegative |
| $T$ | Total modeled duration, $T_n+T_s$ | Minutes; positive |
| $p$ | Fraction of total time spent under stress | 0 to 1 |
| $\bar\mu_i(p)$ | Time-weighted growth rate of strain $i$ | min⁻¹ |
| $a$, $b$ | Mutant-minus-WT growth-rate differences under normal and stressed conditions | min⁻¹ |
| $s(p)$ | Mutant-minus-WT difference in time-weighted growth rate | min⁻¹ |
| $R(t)$ | Mutant-to-WT population ratio, $N_{\mathrm{mut}}(t)/N_{\mathrm{WT}}(t)$ | Dimensionless |
| $p^*$ | Stress time fraction at which $s(p)=0$ | Dimensionless |

The supplied maximum growth rates are treated as constant exponential-growth parameters within each condition. Both strains experience the same environmental schedule. Switching changes the rate immediately, with no lag or memory of the previous condition. Resource limitation, density dependence, and interactions between populations are not modeled. Input rates must be finite and nonnegative, as required by the calculation script.

These assumptions define a comparison of growth potential. The calculation does not estimate survival probability, fixation probability, or uncertainty from biological replicates.

## 2. From population growth to time-weighted growth

For either strain, the starting equation is:

$$
\frac{dN_i(t)}{dt}=\mu_i(t)N_i(t),
\qquad
\mu_i(t)=
\begin{cases}
\mu_{i,n}, & \text{under normal conditions},\\
\mu_{i,s}, & \text{under the selected stress}.
\end{cases}
$$

Dividing by the positive population size gives the instantaneous log-growth rate:

$$
\frac{1}{N_i(t)}\frac{dN_i(t)}{dt}
=\frac{d\ln N_i(t)}{dt}
=\mu_i(t).
$$

Integrating over the full duration gives:

$$
\ln\left(\frac{N_i(T)}{N_i(0)}\right)
=\int_0^T\mu_i(t)\,dt
=\mu_{i,n}T_n+\mu_{i,s}T_s.
$$

Here, $T_n$ and $T_s$ are accumulated exposure times, even if the environment switches repeatedly. Define:

$$
T=T_n+T_s,\qquad
p=\frac{T_s}{T},\qquad
T_n=(1-p)T,\qquad T_s=pT.
$$

Substitution and division by $T$ yield the time-weighted rate:

$$
\bar\mu_i(p)
=\frac{1}{T}\ln\left(\frac{N_i(T)}{N_i(0)}\right)
=(1-p)\mu_{i,n}+p\mu_{i,s}.
$$

Equivalently, the population solution is:

$$
N_i(T)=N_i(0)\exp\left[T\bar\mu_i(p)\right].
$$

The logarithm applies to population fold growth; the measured growth rates themselves are not logarithmically transformed. Under the stated assumptions, only the accumulated exposure times matter, so the order and number of switches do not affect the result at fixed $T_n$ and $T_s$.

If the environment is described by a stress-to-normal duration ratio $r=T_s/T_n$ with $T_n>0$, convert it to the fraction used by the code:

$$
p=\frac{r}{1+r},\qquad r=\frac{p}{1-p}\quad (p<1).
$$

For example, equal stress and normal durations give $r=1$ and $p=0.5$. The endpoint $p=1$ means all time is spent under stress; the duration ratio then has no finite value. The code uses $p$ directly.

## 3. From strain growth rates to relative advantage

First calculate the differences between mutant and WT in each condition:

$$
a=\mu_{\mathrm{mut},n}-\mu_{\mathrm{WT},n},
\qquad
b=\mu_{\mathrm{mut},s}-\mu_{\mathrm{WT},s}.
$$

Subtract their time-weighted rates:

$$
\begin{aligned}
s(p)
&=\bar\mu_{\mathrm{mut}}(p)-\bar\mu_{\mathrm{WT}}(p)\\
&=(1-p)\mu_{\mathrm{mut},n}+p\mu_{\mathrm{mut},s}
-\left[(1-p)\mu_{\mathrm{WT},n}+p\mu_{\mathrm{WT},s}\right]\\
&=(1-p)\left(\mu_{\mathrm{mut},n}-\mu_{\mathrm{WT},n}\right)
+p\left(\mu_{\mathrm{mut},s}-\mu_{\mathrm{WT},s}\right)\\
&=(1-p)a+pb\\
&=a+p(b-a).
\end{aligned}
$$

Thus $s(0)=a$ and $s(1)=b$. Positive $s(p)$ favors the mutant, negative $s(p)$ favors WT, and zero means equal time-weighted rates. The result is a difference in min⁻¹; it is not a growth-rate ratio or a percentage change.

To connect this quantity to population growth, divide the two population solutions:

$$
\begin{aligned}
R(T)
&=\frac{N_{\mathrm{mut}}(T)}{N_{\mathrm{WT}}(T)}\\
&=\frac{N_{\mathrm{mut}}(0)}{N_{\mathrm{WT}}(0)}
\exp\left[T\left(\bar\mu_{\mathrm{mut}}(p)-\bar\mu_{\mathrm{WT}}(p)\right)\right]\\
&=R(0)\exp[T s(p)].
\end{aligned}
$$

Taking the natural logarithm gives:

$$
\ln\left(\frac{R(T)}{R(0)}\right)=T s(p),
\qquad
s(p)=\frac{1}{T}\ln\left(\frac{R(T)}{R(0)}\right).
$$

This is why $s(p)$ is called a relative log-growth rate. At fixed growth rates and fixed $p$, neither $s(p)$ nor the crossover fraction depends on total duration. The accumulated log-ratio change $T s(p)$ and the fold change $\exp[T s(p)]$ do depend on duration. If extending an environmental schedule changes its stress fraction, $p$ must be recalculated.

Equal starting populations give $R(0)=1$. For unequal starting populations, positive $s(p)$ means that the mutant-to-WT ratio increases; it does not necessarily mean the mutant population is larger at time $T$.

## 4. Crossover fraction and advantage direction

An interior crossover exists when the differences have opposite signs: $ab<0$. Solve $s(p^*)=0$:

$$
\begin{aligned}
0&=(1-p^*)a+p^*b\\
 &=a+p^*(b-a),\\
p^*(b-a)&=-a,\\
p^*&=\frac{-a}{b-a}=\frac{a}{a-b}.
\end{aligned}
$$

Opposite signs ensure $b-a\ne0$ and $0<p^*<1$. To determine which side favors the mutant, rewrite the rate difference as:

$$
s(p)=(b-a)(p-p^*).
$$

If $a<0$ and $b>0$, the slope $b-a$ is positive, so the mutant advantage occurs for $p>p^*$. If $a>0$ and $b<0$, the slope is negative, so the advantage occurs for $p<p^*$. At the threshold, the population ratio remains equal to its initial value under the model.

When $a$ and $b$ have the same strict sign, their weighted average has that sign for every $p$ in the interval. The mutant is therefore beneficial throughout the range if both are positive, or deleterious throughout if both are negative. The script leaves `p_star` blank for these cases. If $a=b\ne0$, then $s(p)=a$ is constant and the crossover formula is not used.

The equations also define neutral boundary cases, but the current script rejects them rather than assigning an output category:

| Case | Mathematical result |
| --- | --- |
| $a=0$, $b\ne0$ | $s(p)=pb$: neutral at $p=0$, with the sign of $b$ for $p>0$. |
| $b=0$, $a\ne0$ | $s(p)=(1-p)a$: neutral at $p=1$, with the sign of $a$ for $p<1$. |
| $a=b=0$ | $s(p)=0$ throughout the range. |

The fraction $p^*$ is a crossover in environmental exposure, not a statistical significance threshold.

## 5. Calculation performed by the scripts

The calculation is analytical. The scripts do not numerically integrate population trajectories or select a total simulation duration.

1. Read the 172 input growth rates: WT and 42 mutants under four conditions.
2. For each mutant and each of the three stresses, calculate $a$ and $b$ from the measured rates.
3. Evaluate $s(p)=a+p(b-a)$ at $p=0,0.01,\ldots,1$. This gives 101 scenarios per comparison and 12,726 scenario rows overall.
4. For opposite-sign differences, calculate $p^*$ directly from the formula, independently of the scenario grid. Otherwise, apply the same-sign classification. Endpoint ties stop the calculation with an error.
5. Export the 126 comparison summaries, scenario values, and classification counts. The plotting script uses the same calculation function to obtain its values.

The summary rates at three representative stress fractions are:

$$
s(0.25)=0.75a+0.25b,\qquad
s(0.50)=0.50a+0.50b,\qquad
s(0.75)=0.25a+0.75b.
$$

See the [output reference](README.md#output-files-and-fields) for CSV field names and the [heatmap guide](README.md#reading-the-heatmap) for display units and symbols.

## 6. Worked example with artificial values

Suppose WT has growth rates of 0.010 min⁻¹ under normal conditions and 0.004 min⁻¹ under stress. A mutant has rates of 0.009 and 0.006 min⁻¹, respectively. Then:

$$
a=0.009-0.010=-0.001\ \mathrm{min}^{-1},\qquad
b=0.006-0.004=0.002\ \mathrm{min}^{-1}.
$$

At $p=0.5$, the time-weighted rates are:

$$
\bar\mu_{\mathrm{WT}}=0.5(0.010)+0.5(0.004)=0.007\ \mathrm{min}^{-1},
$$

$$
\bar\mu_{\mathrm{mut}}=0.5(0.009)+0.5(0.006)=0.0075\ \mathrm{min}^{-1}.
$$

Their difference and crossover fraction are:

$$
s(0.5)=0.0075-0.007=0.0005\ \mathrm{min}^{-1},\qquad
p^*=\frac{0.001}{0.003}=\frac{1}{3}.
$$

The mutant has a modeled advantage above a stress time fraction of one third. With equal starting populations and $T=1000$ minutes:

$$
\frac{N_{\mathrm{mut}}(T)}{N_{\mathrm{WT}}(T)}
=\exp(1000\times0.0005)
=\exp(0.5)
\approx1.65.
$$

At $T=2000$ minutes and the same $p$, $s(p)$ remains 0.0005 min⁻¹, while the ratio becomes $\exp(1)\approx2.72$. These artificial examples illustrate the equations and are not experimental results.
