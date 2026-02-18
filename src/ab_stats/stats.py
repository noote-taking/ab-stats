"""
A/B 테스트 통계 검정: 비율 차이(z-test), 평균 차이(Welch t-test)

리스트/관측수를 받아 함수 내부에서 평균·분산·표본크기를 계산한 뒤 검정합니다
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t, norm


def _to_valid_arrays(control_values, treatment_values):
    """
    리스트를 numpy 배열로 변환하고, NaN을 제거
    """
    x = np.asarray(control_values, dtype=float)
    y = np.asarray(treatment_values, dtype=float)
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    return x, y


def proportions_ztest(
    control_n,
    control_success,
    treatment_n,
    treatment_success,
    alpha=0.05,
    power=0.8,
):
    """
    두 독립 비율의 차이에 대한 z-검정 (대조군 vs 실험군).

    대조군/실험군의 관측수와 성공수를 받아, 실험군 비율이 대조군과 다른지 검정하고
    증감률·신뢰구간·최소 샘플 수 등을 반환합니다.

    Parameters
    ----------
    control_n : int
        대조군 관측수.
    control_success : int
        대조군 성공수.
    treatment_n : int
        실험군 관측수.
    treatment_success : int
        실험군 성공수.
    alpha : float, optional
        유의수준 (default: 0.05).
    power : float, optional
        검정력 1 - beta (default: 0.8).

    Returns
    -------
    pandas.DataFrame
        한 행에 metric_formula, metric_value, delta_relative, delta_absolute,
        p_value, CI_relative, CI_absolute, MSS, statistic 컬럼.
    """
    epsilon = 1e-10

    if control_n <= 0 or treatment_n <= 0:
        raise ValueError("control_n and treatment_n must be positive.")
    if not (0 <= control_success <= control_n and 0 <= treatment_success <= treatment_n):
        raise ValueError("Success counts must be between 0 and respective n.")

    p1 = control_success / control_n
    p2 = treatment_success / treatment_n
    diff = p2 - p1

    # 표준오차 및 z-통계량·p-value 직접 계산
    # delta_se: 비율 차이(p2-p1)의 표준오차 (standard error)
    delta_se = np.sqrt(p1 * (1 - p1) / control_n + p2 * (1 - p2) / treatment_n)
    if delta_se <= 0:
        delta_se = epsilon
    stat_z = diff / delta_se  # stat_z: z-통계량
    p_value = 2 * (1 - norm.cdf(abs(stat_z)))  # 양측검정 p-value

    # metric_formula: 실험군 성공수/실험군 관측수
    metric_formula = f"{treatment_success}/{treatment_n}"
    metric_value = p2
    delta_absolute = diff
    delta_relative = (diff / p1) * 100 if p1 > epsilon else np.nan

    # 신뢰구간 (비율 차이: z 사용)
    # z_crit: 신뢰구간 계산을 위한 z-분포의 임계값 (양측검정 상위 α/2 지점)
    z_crit = norm.ppf(1 - alpha / 2)
    delta_ci_lower = diff - z_crit * delta_se  # CI: confidence interval (신뢰구간)
    delta_ci_upper = diff + z_crit * delta_se
    CI_absolute = f"[{delta_ci_lower:.4f}, {delta_ci_upper:.4f}]"

    # 증감률(%) 신뢰구간
    if p1 > epsilon:
        uplift_se = np.sqrt(
            (1 / p1**2) * p2 * (1 - p2) / treatment_n
            + (p2**2 / p1**4) * p1 * (1 - p1) / control_n
        )
        uplift = (diff / p1) * 100
        uplift_ci_lower = uplift - z_crit * uplift_se
        uplift_ci_upper = uplift + z_crit * uplift_se
        CI_relative = f"[{uplift_ci_lower:.2f}%, {uplift_ci_upper:.2f}%]"
    else:
        CI_relative = "[nan%, nan%]"

    # MSS: Minimum Sample Size (최소 샘플 수) 대비 현재 비율
    beta = 1 - power  # beta: 제2종 오류 확률
    z_alpha = norm.ppf(1 - alpha / 2)  # z_alpha: 유의수준 α에 대한 z 임계값
    z_beta = norm.ppf(1 - beta)  # z_beta: 검정력(1-β)에 대한 z 값
    k = control_n / treatment_n if treatment_n > 0 else 0  # k: 대조군/실험군 샘플 비율

    if k <= 0 or not (0 < p1 < 1 and 0 < p2 < 1):
        MSS = "0.0% (∞)"
    else:
        delta_abs = abs(p2 - p1)  # delta_abs: 비율 차이의 절댓값
        if delta_abs <= 0:
            MSS = "0.0% (∞)"
        else:
            # n2_min: 검정력을 만족하기 위한 실험군 최소 샘플 수
            n2_min = ((z_alpha + z_beta) ** 2) * (
                (p1 * (1 - p1)) / k + p2 * (1 - p2)
            ) / (delta_abs**2)
            if np.isnan(n2_min) or n2_min <= 0:
                MSS = "0.0% (∞)"
            else:
                min_n2 = int(np.ceil(n2_min))  # min_n2: 최소 샘플 수 (정수로 올림)
                ratio_pct = (treatment_n / min_n2) * 100  # ratio_pct: 현재 샘플 수 / 최소 샘플 수 (%)
                MSS = f"{ratio_pct:.1f}% ({min_n2:,})"

    return pd.DataFrame(
        [
            {
                "metric_formula": metric_formula,
                "metric_value": metric_value,
                "delta_relative": delta_relative,
                "delta_absolute": delta_absolute,
                "p_value": round(p_value, 5),
                "CI_relative": CI_relative,
                "CI_absolute": CI_absolute,
                "MSS": MSS,
                "statistic": round(float(stat_z), 2),
            }
        ]
    )


def ttest_ind_welch(
    control_values,
    treatment_values,
    alpha=0.05,
    power=0.8,
):
    """
    두 독립 표본의 평균 차이에 대한 Welch t-검정.

    대조군/실험군 **값의 리스트**를 받아, 함수 내부에서 평균·분산·표본크기를 계산한 뒤
    Welch-Satterthwaite 자유도로 검정합니다.

    Parameters
    ----------
    control_values : array_like
        대조군 관측값 리스트 (또는 배열).
    treatment_values : array_like
        실험군 관측값 리스트 (또는 배열).
    alpha : float, optional
        유의수준 (default: 0.05).
    power : float, optional
        검정력 1 - beta (default: 0.8).

    Returns
    -------
    pandas.DataFrame
        한 행에 metric_formula, metric_value, delta_relative, delta_absolute,
        p_value, CI_relative, CI_absolute, MSS, statistic, df 컬럼.
    """
    x, y = _to_valid_arrays(control_values, treatment_values)
    n1, n2 = len(x), len(y)

    if n1 <= 1 or n2 <= 1:
        raise ValueError(
            "Each group must have at least 2 observations (for variance). "
            f"Got control n={n1}, treatment n={n2}."
        )

    # 리스트 기준: 관측수 = len, 누적값 = sum (0 포함)
    treatment_sum = float(np.sum(y))
    mu1 = float(np.mean(x))
    mu2 = float(np.mean(y))
    # 표본 분산 (ddof=1)
    s1_sq = float(np.var(x, ddof=1)) if n1 > 1 else 0.0
    s2_sq = float(np.var(y, ddof=1)) if n2 > 1 else 0.0

    # metric_formula: 실험군 누적값/실험군 관측수 (분자는 정수로 표기, 합계가 소수인 경우 없음)
    metric_formula = f"{int(treatment_sum)}/{n2}"
    metric_value = mu2

    # Welch t-통계량 및 자유도
    se = np.sqrt(s1_sq / n1 + s2_sq / n2)  # se: standard error (표준오차)
    if se <= 0:
        se = 1e-10
    t_stat = (mu2 - mu1) / se  # t_stat: t-통계량

    # Welch-Satterthwaite 자유도 계산: 분자(num_df)와 분모(den_df)
    num_df = (s1_sq / n1 + s2_sq / n2) ** 2  # num_df: numerator of degrees of freedom
    den_df = (s1_sq / n1) ** 2 / (n1 - 1) + (s2_sq / n2) ** 2 / (n2 - 1)  # den_df: denominator
    df_welch = num_df / den_df if den_df > 0 else np.nan  # df_welch: Welch 자유도

    # 두 그룹 모두 분산 0(상수)이면 자유도·p-value·CI 무의미 → NaN 반환
    if np.isfinite(df_welch) and df_welch > 0:
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_welch))
        # t_crit_for_ci: 신뢰구간(CI) 계산을 위한 t-분포의 임계값 (critical value)
        # 양측검정에서 상위 α/2 지점의 t 값 (예: α=0.05면 97.5% 지점)
        t_crit_for_ci = t.ppf(1 - alpha / 2, df_welch)
    else:
        p_value = np.nan
        t_crit_for_ci = np.nan

    diff = mu2 - mu1
    delta_absolute = diff
    epsilon = 1e-10
    delta_relative = (diff / mu1) * 100 if abs(mu1) > epsilon else np.nan

    # 신뢰구간 (평균 차이: t 사용)
    delta_ci_lower = diff - t_crit_for_ci * se
    delta_ci_upper = diff + t_crit_for_ci * se
    CI_absolute = f"[{delta_ci_lower:.4f}, {delta_ci_upper:.4f}]"

    # 증감률(%) 신뢰구간
    if abs(mu1) > epsilon:
        # uplift_se: 증감률(uplift = (μ2-μ1)/μ1 * 100%)의 표준오차
        uplift_se = np.sqrt(
            (1 / max(mu1, epsilon) ** 2) * s2_sq / n2
            + (mu2**2 / max(mu1, epsilon) ** 4) * s1_sq / n1
        )
        uplift = (diff / mu1) * 100  # uplift: 증감률 (%)
        uplift_ci_lower = uplift - t_crit_for_ci * uplift_se
        uplift_ci_upper = uplift + t_crit_for_ci * uplift_se
        CI_relative = f"[{uplift_ci_lower:.2f}%, {uplift_ci_upper:.2f}%]"
    else:
        CI_relative = "[nan%, nan%]"

    # MSS: Minimum Sample Size (최소 샘플 수) - 평균 차이 기준 (Welch 가정, k = n1/n2)
    beta = 1 - power  # beta: 제2종 오류 확률
    z_alpha = norm.ppf(1 - alpha / 2)  # z_alpha: 유의수준 α에 대한 z 임계값
    z_beta = norm.ppf(1 - beta)  # z_beta: 검정력(1-β)에 대한 z 값
    k = n1 / n2 if n2 > 0 else 0  # k: 대조군/실험군 샘플 비율

    if k <= 0:
        MSS = "0.0% (∞)"
    else:
        delta_abs = abs(mu2 - mu1)  # delta_abs: 평균 차이의 절댓값
        if delta_abs <= 0 or np.isnan(delta_abs):
            MSS = "0.0% (∞)"
        else:
            # n2_min: 검정력을 만족하기 위한 실험군 최소 샘플 수
            n2_min = ((z_alpha + z_beta) ** 2) * (s1_sq / k + s2_sq) / (delta_abs**2)
            if np.isnan(n2_min) or n2_min <= 0:
                MSS = "0.0% (∞)"
            else:
                min_n2 = int(np.ceil(n2_min))  # min_n2: 최소 샘플 수 (정수로 올림)
                ratio_pct = (n2 / min_n2) * 100  # ratio_pct: 현재 샘플 수 / 최소 샘플 수 (%)
                MSS = f"{ratio_pct:.1f}% ({min_n2:,})"

    return pd.DataFrame(
        [
            {
                "metric_formula": metric_formula,
                "metric_value": metric_value,
                "delta_relative": delta_relative,
                "delta_absolute": delta_absolute,
                "p_value": round(p_value, 5),
                "CI_relative": CI_relative,
                "CI_absolute": CI_absolute,
                "MSS": MSS,
                "statistic": round(float(t_stat), 2),
                "df": round(float(df_welch), 2),
            }
        ]
    )
