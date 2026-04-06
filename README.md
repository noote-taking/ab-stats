# ab-stats

**ab-stats** is a lightweight Python library for **A/B test statistical analysis**.
It provides a two-sample proportion z-test and Welch's t-test with confidence intervals,
uplift, and post-hoc minimum sample size in a pandas DataFrame. For **continuous metrics** with several experiment groups, you can **optionally** use **`winsorize_experiment_groups`** first to pool all groups, compute percentiles, and apply upper-tail caps before pairwise Welch tests.


## When to use which test?

Choose the function that matches your metric: a **proportion (rate)** or a **mean (continuous value)** to compare between control and treatment groups.

| Metric type       | Example                                     | Function              |
|-------------------|---------------------------------------------|-----------------------|
| Proportion (rate) | CTR, PUR, conversion rate, signup rate      | `proportions_ztest()` |
| Mean (continuous) | ARPU, average session length, time on page  | `ttest_ind_welch()`   |

## Key notes
- **Rich output**: `proportions_ztest` and `ttest_ind_welch` each return a one-row DataFrame with `control_formula`, `treatment_formula`, `control_value`, `treatment_value`, `delta_relative`, `delta_absolute`, `p_value`, `CI_relative`, `CI_absolute`, `MSS_posthoc`, and `statistic`. **`ttest_ind_welch` also includes `df`.**
- **Two-sided tests**: Both tests perform two-sided hypothesis tests
- **Delta method**: Confidence intervals for uplift (relative change) are computed using the delta method
- **Note on MSS_posthoc**: Minimum sample size required for the given α and β under the assumption that the observed effect is true. **It is computed post hoc and is for reference only** (not for pre-experiment sample size calculation).


## Installation

### Dependencies

ab-stats depends on:

- NumPy (>= 1.20)
- Pandas (>= 1.3)
- SciPy (>= 1.7)

Python 3.8 or newer is required.

### User installation

Install from PyPI with pip:

```bash
pip install ab-stats
```

or with conda (Conda packages are planned):

```bash
conda install -c conda-forge ab-stats
```

## Quick start

### 1. proportions_ztest()

Use this when your metric is a **rate** (e.g. conversion rate, click-through rate). Pass **sample sizes** and **success counts** for control and treatment; the function returns uplift, confidence intervals, and minimum sample size.

> **Parameters**

- **`control_n`** : *int*
  Total number of observations in the control group.

- **`control_success`** : *int*  
  Number of “successes” (e.g. converted users) in the control group.

- **`treatment_n`** : *int*  
  Total number of observations in the treatment group.

- **`treatment_success`** : *int*
  Number of “successes” in the treatment group.

- **`alpha`** : *float, optional*  
  Significance level for confidence intervals and MSS_posthoc.  
  Default is `0.05` (95% confidence interval).

- **`power`** : *float, optional*  
  Target statistical power 1 − β used when computing MSS_posthoc.  
  Default is `0.8` (80% power).

> **Returns**  
> *pandas.DataFrame (one row) with the following columns:*

  - **`control_formula`** : *str*  
    Control metric formula (e.g. `101/998`).

  - **`treatment_formula`** : *str*  
    Treatment metric formula (e.g. `122/1001`).

  - **`control_value`** : *str*  
    Control proportion, formatted as percentage (e.g. `10.12%`).

  - **`treatment_value`** : *str*  
    Treatment proportion, formatted as percentage (e.g. `12.19%`).

  - **`delta_relative`** : *str*  
    Relative change (uplift) of treatment vs control, formatted as a percentage (e.g. `20.43%`).

  - **`delta_absolute`** : *str*  
    Absolute difference in proportions in percentage points (e.g. `2.07%p`).

  - **`p_value`** : *float*  
    Two-sided p-value for the null hypothesis p₁ = p₂.

  - **`CI_relative`** : *str*  
    Confidence interval for the relative change (uplift), formatted as `[L%, U%]`.

  - **`CI_absolute`** : *str*  
    Confidence interval for the absolute difference in proportions, formatted as `[L, U]`.

  - **`MSS_posthoc`** : *str*  
    Post hoc minimum sample size status (e.g. `27.5% (3,641)`).  
    The percentage is the ratio of current treatment sample size to the required minimum under the observed effect.

  - **`statistic`** : *float*  
    z-statistic of the two-sample proportion test.


#### Example
```python
from ab_stats import proportions_ztest

# Control: 101 successes out of 998; Treatment: 122 successes out of 1001
df = proportions_ztest(
    control_n=998,
    control_success=101,
    treatment_n=1001,
    treatment_success=122,
    alpha=0.05,
    power=0.8,
)
print(df)
```

**Output:**

| control_formula | treatment_formula | control_value | treatment_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS_posthoc | statistic |
|-----------------|-------------------|---------------|-----------------|----------------|----------------|---------|-------------|-------------|-------------|-----------|
| 101/998 | 122/1001 | 10.12% | 12.19% | 20.43% | 2.07%p | 0.14206 | [-9.52%, 50.38%] | [-0.01, 0.05] | 27.5% (3,641) | 1.47 |

### 2. ttest_ind_welch()

Use this when your metric is a **mean** (e.g. average revenue per user, average session length). Pass **lists of values** (one value per user or per observation) for control and treatment; the function computes means, variances, and sample sizes internally and returns uplift, confidence intervals, and minimum sample size. The result includes *df* (Welch–Satterthwaite degrees of freedom). **Winsorization is not applied inside this function**; use `winsorize_experiment_groups` (section 3 below) first if you need pool-based upper winsorization before Welch.

> **Parameters**

- **`control_values`** : *array_like*  
  Observations in the control group (e.g. list or array; one value per observation).

- **`treatment_values`** : *array_like*  
  Observations in the treatment group.

- **`alpha`** : *float, optional*
  Significance level for confidence intervals and MSS_posthoc.  
  Default is `0.05` (95% confidence interval).

- **`power`** : *float, optional*
  Target statistical power 1 − β used when computing MSS_posthoc.  
  Default is `0.8` (80% power).

> **Returns**  
> *pandas.DataFrame (one row) with the following columns:*

  - **`control_formula`** : *str*  
    Control mean formula: sum of observations / *n*, with the sum shown to two decimal places (e.g. `123.60/12`).

  - **`treatment_formula`** : *str*  
    Treatment mean formula (e.g. `107.80/10`).

  - **`control_value`** : *float*  
    Observed mean in the control group.

  - **`treatment_value`** : *float*  
    Observed mean in the treatment group.

  - **`delta_relative`** : *str*  
    Relative change (uplift) of treatment vs control mean, formatted as a percentage.

  - **`delta_absolute`** : *float*  
    Absolute difference in means (treatment − control).

  - **`p_value`** : *float*  
    Two-sided p-value for the null hypothesis μ₁ = μ₂.

  - **`CI_relative`** : *str*  
    Confidence interval for the relative change (uplift), formatted as `[L%, U%]`.

  - **`CI_absolute`** : *str*  
    Confidence interval for the absolute difference in means, formatted as `[L, U]`.

  - **`MSS_posthoc`** : *str*  
    Post hoc minimum sample size status (e.g. `62.5% (16)`).  
    The percentage is the ratio of current treatment sample size to the required minimum under the observed effect.

  - **`statistic`** : *float*  
    t-statistic of Welch’s t-test.

  - **`df`** : *float*  
    Degrees of freedom (Welch–Satterthwaite).

#### Example (raw values, no winsorization)
```python
from ab_stats import ttest_ind_welch

control = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment = [11.0, 10.5, 11.8, 10.9, 11.2, 10.5, 10.7, 10.1, 10.3, 10.8]  # n=10

df = ttest_ind_welch(control_values=control, treatment_values=treatment, alpha=0.05, power=0.8)
print(df)
```

**Output:**

| control_formula | treatment_formula | control_value | treatment_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS_posthoc | statistic | df |
|-----------------|-------------------|---------------|-----------------|----------------|----------------|---------|-------------|-------------|-------------|-----------|-----|
| 123.60/12 | 107.80/10 | 10.3 | 10.78 | 4.66% | 0.48 | 0.03383 | [0.30%, 9.02%] | [0.04, 0.92] | 62.5% (16) | 2.28 | 19.41 |

### 3. winsorize_experiment_groups()

**Optional.** For **continuous** metrics, call this **before** `ttest_ind_welch` when you want pool-based **upper** winsorization. It pools **all** groups into one array, computes pool percentiles (and an optional cap from the pool), then applies that cap **per group**—one threshold for the whole experiment.

> **Parameters**

- **`observations_by_variant`** : *mapping of str to array_like*  
  **Keys** are variant labels (e.g. `"A"`, `"B"`). **Values** are 1-D numeric arrays: all metric values for that variant (typically one value per user).  
  Build this from a long-format SQL table by grouping on `variant_name` and collecting the metric column (e.g. `cum_event_value`) into `{"A": ndarray, "B": ndarray, ...}`.

- **`upper_tail`** : *{None, 0.05, 0.01}, optional*  
  Which **upper-tail** rule sets the pool-based cap: `None` copies each group without capping; `0.05` caps at the pool **95th** percentile; `0.01` at the pool **99th** percentile.

> **Returns**  
> *`WinsorizeExperimentResult`* with fields:

  - **`percentiles`** : *str* — formatted 90th–100th pool percentiles (100th = pool maximum).  
  - **`p90`, `p95`, `p99`, `p100`** : *float* — same numbers as numeric values.  
  - **`pool_n`** : *int* — total count of non-NaN values across all experiment groups.  
  - **`winsorized`** : *dict[str, numpy.ndarray]* — same keys as `observations_by_variant`, capped arrays.

#### Example 1 — Two experiment groups

**1. Winsorize (pool percentiles, then cap per group)**

```python
from ab_stats import winsorize_experiment_groups, ttest_ind_welch

control = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment = [11.0, 10.5, 11.8, 10.9, 11.2, 10.5, 10.7, 10.1, 10.3, 10.8]  # n=10

winsorize_result = winsorize_experiment_groups(
    {"A": control, "B": treatment}, upper_tail=0.01
)
print(winsorize_result.percentiles)
```

**Output:**

```
[90th: 11.18, 95th: 11.20, 99th: 11.67, 100th: 11.80]
```

**2. Welch on `winsorize_result.winsorized["A"]` vs `["B"]`**

```python
df_ab = ttest_ind_welch(
    control_values = winsorize_result.winsorized["A"],
    treatment_values = winsorize_result.winsorized["B"],
    alpha=0.05,
    power=0.8,
)
print(df_ab)
```

**Output:**

| control_formula | treatment_formula | control_value | treatment_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS_posthoc | statistic | df |
|-----------------|-------------------|---------------|-----------------|----------------|----------------|---------|-------------|-------------|-------------|-----------|-----|
| 123.60/12 | 107.67/10 | 10.3 | 10.7674 | 4.54% | 0.467 | 0.03286 | [0.32%, 8.76%] | [0.04, 0.89] | 66.7% (15) | 2.29 | 19.74 |

#### Example 2 — Three or more experiment groups

In production, **`exp_values_by_variant`** should be the **long-format** query table as a `pandas.DataFrame` (e.g. columns `variant_name`, `cum_event_value`). The snippet below uses that name for a small stand-in `DataFrame` so the **Output** block matches when you run it.

```python
import pandas as pd
from ab_stats import winsorize_experiment_groups, ttest_ind_welch

# Arbitrary toy rows: long-format table with variant_name + one metric column (shape you need from SQL).
exp_values_by_variant = pd.DataFrame(
    {
        "variant_name": ["A", "A", "A", "B", "B", "B", "C", "C", "C"],
        "cum_event_value": [10.1, 9.8, 11.2, 11.0, 10.5, 11.8, 10.0, 10.2, 10.4],
    }
)

variant_names = sorted(exp_values_by_variant["variant_name"].dropna().unique().tolist())
print("variant_names:", variant_names)

event_values_by_variant = {
    variant_name: exp_values_by_variant.loc[
        exp_values_by_variant["variant_name"] == variant_name, "cum_event_value"
    ].astype(float).values
    for variant_name in variant_names
}

winsorize_result = winsorize_experiment_groups(event_values_by_variant, upper_tail=0.01)
print("percentiles:", winsorize_result.percentiles)
print("pool_n:", winsorize_result.pool_n)

df_ab = ttest_ind_welch(
    control_values=winsorize_result.winsorized["A"],
    treatment_values=winsorize_result.winsorized["B"],
    alpha=0.05,
    power=0.8,
)
print(df_ab.to_string(index=False))
```

Pairwise comparisons for other groups (e.g. `A` vs `C`) use the same pattern with the corresponding keys on `winsorize_result.winsorized`.

**Output:**

```
variant_names: ['A', 'B', 'C']
percentiles: [90th: 11.32, 95th: 11.56, 99th: 11.75, 100th: 11.80]
pool_n: 9
control_formula treatment_formula  control_value  treatment_value delta_relative  delta_absolute  p_value      CI_relative   CI_absolute MSS_posthoc  statistic   df
        31.10/3           33.25/3      10.366667           11.084          6.92%           0.717  0.27092 [-8.83%, 22.67%] [-0.85, 2.29]  20.0% (15)       1.28 3.91
```

### 4. Using with Pandas

Results are returned as a pandas DataFrame, so you can merge with other columns or filter as usual.

```python
from ab_stats import proportions_ztest, ttest_ind_welch

# Proportion test
result_prop = proportions_ztest(1000, 100, 1000, 120)
print("Proportion test:")
print("MSS_posthoc: ", result_prop["MSS_posthoc"].iloc[0])
print("treatment_value: ", result_prop["treatment_value"].iloc[0])
print("delta_relative: ", result_prop["delta_relative"].iloc[0])
print("p_value: ", result_prop["p_value"].iloc[0])
print("CI_relative: ", result_prop["CI_relative"].iloc[0])
print("statistic: ", result_prop["statistic"].iloc[0])

# Mean test
control_vals = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment_vals = [11.0, 11.5, 11.8, 11.9, 11.2, 11.5, 10.7, 11.1, 10.3, 10.8]  # n=10
result_ttest = ttest_ind_welch(control_vals, treatment_vals)
print("\nMean test:")
print("control_value: ", result_ttest["control_value"].iloc[0])
print("delta_relative: ", result_ttest["delta_relative"].iloc[0])
print("p_value: ", result_ttest["p_value"].iloc[0])
print("df: ", result_ttest["df"].iloc[0])
```

**Output:**

```
Proportion test:
MSS_posthoc:  26.0% (3,839)
treatment_value:  12.00%
delta_relative:  20.00%
p_value:  0.15292
CI_relative:  [-10.06%, 50.06%]
statistic:  1.43

Mean test:
control_value:  10.3
delta_relative:  8.54%
p_value:  0.0006
df:  19.15
```

## References

- [1] Zhou, J., Lu, J., & Shallah, A. (2023). All about sample-size calculations for A/B testing: Novel extensions & practical guide. *Proceedings of the 32nd ACM International Conference on Information and Knowledge Management (CIKM '23)*, 1–30.
- [2] Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). *Sample Size Calculations in Clinical Research* (3rd ed.). Chapman & Hall/CRC Biostatistics Series.
- [3] noote-taking. (n.d.). When and how to calculate minimum sample size. noote-taking.github.io. https://noote-taking.github.io/%ED%86%B5%EA%B3%84%ED%95%99/when-and-how-to-calculate-minimum-sample-size/

## License

MIT License. See [LICENSE](LICENSE) for details.
