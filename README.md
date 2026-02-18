# ab-stats

**ab-stats** is a Python library that computes the statistics you need for A/B tests. It runs a **two-sample proportion z-test** for **rate (proportion) differences** and **Welch's t-test** for **mean differences** between control and treatment groups, and returns p-value, confidence intervals, uplift (relative change), and minimum sample size (MSS) in a pandas DataFrame.

## Documentation

- [README](https://github.com/noote-taking/ab-stats#readme)

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

or with conda:

```bash
conda install -c conda-forge ab-stats
```

## Quick start

### 1. Proportion (rate) difference — `proportions_ztest`

Pass **sample sizes** and **success counts** for control and treatment; the function runs a two-sample proportion z-test and returns uplift, confidence intervals, and minimum sample size. MSS is the sample size required for the given α and β under the assumption that the observed effect is true; it is computed post hoc and should be used as a reference only.

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
| 122/1001 | 0.1219 | 20.45 | 0.0207 | 0.14162 | [5.12%, 35.78%] | [-0.0069, 0.0483] | 152.3% (657) | 1.47 |

### 2. Mean difference — `ttest_ind_welch`

Pass **lists of values** for control and treatment; the function computes means, variances, and sample sizes internally and runs Welch's t-test.

```python
from ab_stats import ttest_ind_welch

# Example: observation lists for control and treatment
control = [10.1, 9.8, 11.2, 10.5, 9.9, 10.8, 10.3, 11.0, 9.7, 10.4, 9.8, 10.1]  # n=12
treatment = [12.0, 11.5, 12.8, 11.9, 12.2, 12.5, 11.7, 12.1, 12.3, 11.8]  # n=10

df = ttest_ind_welch(control, treatment, alpha=0.05, power=0.8)
print(df)
```

**Output:**

| metric_formula | metric_value | delta_relative | delta_absolute | p_value | CI_relative | CI_absolute | MSS | statistic | df |
|----------------|--------------|----------------|----------------|---------|-------------|-------------|-----|-----------|-----|
| 120/10 | 12.03 | 17.14 | 1.76 | 0.00273 | [8.21%, 26.07%] | [0.65, 2.87] | 45.2% (221) | 3.45 | 18.52 |

### 3. Using with Pandas

Results are returned as a pandas DataFrame, so you can merge with other columns or filter as usual.

```python
import pandas as pd
from ab_stats import proportions_ztest, ttest_ind_welch

# Proportion test
result_prop = proportions_ztest(1000, 100, 1000, 120)
# Use result_prop["p_value"], result_prop["CI_relative"], etc.

# Mean test (lists → means, variances, n are computed inside the function)
control_vals = [1.0, 2.0, 3.0, 4.0, 5.0]
treatment_vals = [2.0, 3.0, 4.0, 5.0, 6.0]
result_ttest = ttest_ind_welch(control_vals, treatment_vals)
# Use result_ttest["metric_value"], result_ttest["df"], etc.
```

## License

MIT License. See [LICENSE](LICENSE) for details.
