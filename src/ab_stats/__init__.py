"""
ab-stats: A/B test statistical testing library.

Two-sample proportion z-test and Welch t-test with p-value, confidence
intervals, uplift, and minimum sample size (MSS).
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("ab-stats")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

from ab_stats.stats import proportions_ztest, ttest_ind_welch

__all__ = ["proportions_ztest", "ttest_ind_welch", "__version__"]
