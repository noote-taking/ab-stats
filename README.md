# ab-stats

**ab-stats** is a Python library that computes the statistics you need for A/B tests. It runs a **two-sample proportion z-test** for **rate (proportion) differences** and **Welch's t-test** for **mean differences** between control and treatment groups, and returns p-value, confidence intervals, uplift (relative change), and minimum sample size (MSS) in a pandas DataFrame.

## Features
- `proportions_ztest()`: Tests the difference in proportion (rate) metrics between control and treatment groups
- `ttest_ind_welch()`: Tests the difference in mean metrics between control and treatment groups

## Key notes
- **Rich output**: Returns pandas DataFrame with metric_formula, metric_value, delta_relative, delta_absolute, p_value, CI_relative, CI_absolute, MSS, statistic (and df for t-test)
- **Two-sided tests**: Both functions perform two-sided hypothesis tests
- **Delta method**: Confidence intervals for uplift (relative change) computed using the delta method
- **Note on MSS**: Minimum Sample Size (MSS) is the sample size required for the given α and β under the assumption that the observed effect is true. It is computed post hoc and should be used as a reference only (applies to both proportion and mean tests).

## Documentation

- [GitHub Repository](https://github.com/noote-taking/ab-stats)

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

or with conda (currently under review):

```bash
conda install -c conda-forge ab-stats
```

## Quick start

### 1. Proportion (rate) difference — `proportions_ztest()`

Use this when your metric is a **rate** (e.g. conversion rate, click-through rate). Pass **sample sizes** and **success counts** for control and treatment; the function returns uplift, confidence intervals, and minimum sample size.

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

| metric_formula | metric_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS | statistic |
|----------------|--------------|----------------|----------------|---------|-------------|-------------|-----|-----------|
| 122/1001 | 0.121878 | 20.430065 | 0.020676 | 0.1418 | [20.13%, 20.73%] | [-0.0069, 0.0483] | 27.5% (3,641) | 1.47 |

### 2. Mean difference — `ttest_ind_welch()`

Use this when your metric is a **mean** (e.g. average revenue per user, average session length). Pass **lists of values** (one value per user or per observation) for control and treatment; the function computes means, variances, and sample sizes internally and returns uplift, confidence intervals, and minimum sample size. The result also includes *df* (degrees of freedom for the Welch t-test).

```python
from ab_stats import ttest_ind_welch

# Example: observation lists for control and treatment
control = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment = [11.0, 10.5, 11.8, 10.9, 11.2, 10.5, 10.7, 10.1, 10.3, 10.8]  # n=10

df = ttest_ind_welch(control, treatment, alpha=0.05, power=0.8)
print(df)
```

**Output:**

| metric_formula | metric_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS | statistic | df |
|----------------|--------------|----------------|----------------|---------|-------------|-------------|-----|-----------|-----|
| 107/10 | 10.78 | 4.660194 | 0.48 | 0.03383 | [4.62%, 4.70%] | [0.0406, 0.9194] | 62.5% (16) | 2.28 | 19.41 |

### 3. Using with Pandas

Results are returned as a pandas DataFrame, so you can merge with other columns or filter as usual.

```python
from ab_stats import proportions_ztest, ttest_ind_welch

# Proportion test
result_prop = proportions_ztest(1000, 100, 1000, 120)
print("Proportion test:")
print("MSS: ", result_prop["MSS"].iloc[0])
print("metric_value: ", result_prop["metric_value"].iloc[0])
print("delta_relative: ", result_prop["delta_relative"].iloc[0], "%")
print("p_value: ", result_prop["p_value"].iloc[0])
print("CI_relative: ", result_prop["CI_relative"].iloc[0])
print("statistic: ", result_prop["statistic"].iloc[0])

# Mean test
control_vals = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment_vals = [12.0, 11.5, 12.8, 11.9, 12.2, 12.5, 11.7, 12.1, 12.3, 11.8]  # n=10
result_ttest = ttest_ind_welch(control_vals, treatment_vals)
print("\nMean test:")
print("metric_value: ", result_ttest["metric_value"].iloc[0])
print("delta_relative: ", result_ttest["delta_relative"].iloc[0], "%")
print("p_value: ", result_ttest["p_value"].iloc[0])
print("df: ", result_ttest["df"].iloc[0])
```

**Output:**

```
Proportion test:
MSS:  26.0% (3,839)
metric_value:  0.12
delta_relative:  19.99999999999999 %
p_value:  0.15271
CI_relative:  [19.70%, 20.30%]
statistic:  1.43

Mean test:
metric_value:  12.08
delta_relative:  17.281553398058268 %
p_value:  0.0
df:  19.95
```

## License

MIT License. See [LICENSE](LICENSE) for details.
