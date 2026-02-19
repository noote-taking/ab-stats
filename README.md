# ab-stats

**ab-stats** is a lightweight Python library for **A/B test statistical analysis**.
It provides a two-sample proportion z-test and Welch's t-test with confidence intervals,
uplift, and post-hoc minimum sample size in a pandas DataFrame.


## When to use which test?

Use the appropriate function depending on whether your metric is a proportion (rate) or a mean, i.e., choose the test that matches the scale of the metric you want to compare between control and treatment groups.

| Metric type       | Example                                     | Function              |
|-------------------|---------------------------------------------|-----------------------|
| Proportion (rate) | CTR, PUR, conversion rate, signup rate      | `proportions_ztest()` |
| Mean (continuous) | ARPU, average session length, time on page  | `ttest_ind_welch()`   |

## Key notes
- **Rich output**: Returns pandas DataFrame with `metric_formula`, `metric_value`, `delta_relative`, `delta_absolute`,  
  `p_value`, `CI_relative`, `CI_absolute`, `MSS_posthoc`, `statistic` (and `df` for t-test)
- **Two-sided tests**: Both functions perform two-sided hypothesis tests
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
  Target statistical power \(1 - \beta\) used when computing MSS_posthoc.  
  Default is `0.8` (80% power).

> **Returns**  
> *pandas.DataFrame (one row) with the following columns:*

  - **`metric_formula`** : *str*  
    String representation of the metric (e.g. `122/1001`).

  - **`metric_value`** : *float*  
    Observed proportion in the treatment group.

  - **`delta_relative`** : *str*  
    Relative change (uplift) of treatment vs control, formatted as a percentage (e.g. `20.43%`).

  - **`delta_absolute`** : *float*  
    Absolute difference in proportions (treatment − control).

  - **`p_value`** : *float*  
    Two-sided p-value for the null hypothesis \(p_1 = p_2\).

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

| metric_formula | metric_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS_posthoc | statistic |
|----------------|--------------|----------------|----------------|---------|-------------|-------------|-------------|-----------|
| 122/1001 | 0.121878 | 20.43% | 0.02 | 0.1418 | [-9.52%, 50.38%] | [-0.01, 0.05] | 27.5% (3,641) | 1.47 |

### 2. ttest_ind_welch()

Use this when your metric is a **mean** (e.g. average revenue per user, average session length). Pass **lists of values** (one value per user or per observation) for control and treatment; the function computes means, variances, and sample sizes internally and returns uplift, confidence intervals, and minimum sample size. The result also includes *df* (degrees of freedom for the Welch t-test).

> **Parameters**

- **`control_values`** : *array_like*  
  Observations in the control group (e.g. list or array; one value per observation).

- **`treatment_values`** : *array_like*  
  Observations in the treatment group.

- **`alpha`** : *float, optional*
  Significance level for confidence intervals and MSS_posthoc.  
  Default is `0.05` (95% confidence interval).

- **`power`** : *float, optional*
  Target statistical power \(1 - \beta\) used when computing MSS_posthoc.  
  Default is `0.8` (80% power).

> **Returns**  
> *pandas.DataFrame (one row) with the following columns:*

  - **`metric_formula`** : *str*  
    String representation of the treatment mean (e.g. `107/10`, i.e. sum / n).

  - **`metric_value`** : *float*  
    Observed mean in the treatment group.

  - **`delta_relative`** : *str*  
    Relative change (uplift) of treatment vs control mean, formatted as a percentage.

  - **`delta_absolute`** : *float*  
    Absolute difference in means (treatment − control).

  - **`p_value`** : *float*  
    Two-sided p-value for the null hypothesis \(\mu_1 = \mu_2\).

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

#### Example
```python
from ab_stats import ttest_ind_welch

# Example: observation lists for control and treatment
control = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment = [11.0, 10.5, 11.8, 10.9, 11.2, 10.5, 10.7, 10.1, 10.3, 10.8]  # n=10

df = ttest_ind_welch(control_values=control, treatment_values=treatment, alpha=0.05, power=0.8)
print(df)
```

**Output:**

| metric_formula | metric_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS_posthoc | statistic | df |
|----------------|--------------|----------------|----------------|---------|-------------|-------------|-------------|-----------|-----|
| 107/10 | 10.78 | 4.66% | 0.48 | 0.03383 | [0.30%, 9.02%] | [0.04, 0.92] | 62.5% (16) | 2.28 | 19.41 |

### 3. Using with Pandas

Results are returned as a pandas DataFrame, so you can merge with other columns or filter as usual.

```python
from ab_stats import proportions_ztest, ttest_ind_welch

# Proportion test
result_prop = proportions_ztest(1000, 100, 1000, 120)
print("Proportion test:")
print("MSS_posthoc: ", result_prop["MSS_posthoc"].iloc[0])
print("metric_value: ", result_prop["metric_value"].iloc[0])
print("delta_relative: ", result_prop["delta_relative"].iloc[0])
print("p_value: ", result_prop["p_value"].iloc[0])
print("CI_relative: ", result_prop["CI_relative"].iloc[0])
print("statistic: ", result_prop["statistic"].iloc[0])

# Mean test
control_vals = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment_vals = [11.0, 11.5, 11.8, 11.9, 11.2, 11.5, 10.7, 11.1, 10.3, 10.8]  # n=10
result_ttest = ttest_ind_welch(control_vals, treatment_vals)
print("\nMean test:")
print("metric_value: ", result_ttest["metric_value"].iloc[0])
print("delta_relative: ", result_ttest["delta_relative"].iloc[0])
print("p_value: ", result_ttest["p_value"].iloc[0])
print("df: ", result_ttest["df"].iloc[0])
```

**Output:**

```
Proportion test:
MSS_posthoc:  26.0% (3,839)
metric_value:  0.12
delta_relative:  20.00%
p_value:  0.15271
CI_relative:  [-10.06%, 50.06%]
statistic:  1.43

Mean test:
metric_value:  11.18
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
