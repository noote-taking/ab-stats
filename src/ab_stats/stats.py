"""
A/B 테스트 통계 검정: 비율 차이(z-test), 평균 차이(Welch t-test)

리스트/관측수를 받아 함수 내부에서 평균·분산·표본크기를 계산한 뒤 검정합니다
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t, norm


def _strip_nan_float_array(values) -> np.ndarray:
    """array_like → float ndarray, NaN 제거."""
    float_values = np.asarray(values, dtype=float)
    return float_values[~np.isnan(float_values)]


@dataclass(frozen=True)
class WinsorizeExperimentResult:
    """
    여러 실험그룹을 합친 풀 기준 winsorization 결과.

    Attributes
    ----------
    percentiles : str
        합동 풀의 90·95·99·100번째 백분위수 문자열 (100th는 최댓값).
    p90, p95, p99, p100 : float
        동일 풀의 수치 백분위수.
    pool_n : int
        합동 풀 관측 수 (모든 실험그룹의 유효 값 개수).
    winsorized : dict[str, numpy.ndarray]
        입력과 같은 키, 실험그룹별 1차원 배열 (upper_tail=None 이면 원본 복사).
    """

    percentiles: str
    p90: float
    p95: float
    p99: float
    p100: float
    pool_n: int
    winsorized: Dict[str, np.ndarray]


def winsorize_experiment_groups(
    observations_by_variant: Mapping[str, Any],
    upper_tail=None,
) -> WinsorizeExperimentResult:
    """
    모든 실험그룹의 값을 **하나의 합동 풀**로 이어 붙여 백분위를 구합니다.
    `upper_tail`이 0.05 또는 0.01이면 풀에서 정한 **상단 캡**을 각 그룹 배열에 동일하게 적용하고, `None`이면 복사만 합니다.
    풀의 수치에는 키(실험그룹 이름)가 쓰이지 않습니다.

    Parameters
    ----------
    observations_by_variant : mapping of str to array_like
        실험그룹 이름(예: "A", "B")을 키로, 해당 그룹의 지표값들을 1차원 배열로 담은 매핑입니다.
    upper_tail : {None, 0.05, 0.01}, optional
        상단 꼬리 기준으로 풀의 어느 분위수를 캡으로 쓸지 정합니다.
        `None`이면 상단 캡 없이 복사만 합니다.
        `0.05`이면 합동 풀의 **95분위**를 상단 캡으로 사용합니다(상단 약 5% 꼬리).
        `0.01`이면 합동 풀의 **99분위**를 상단 캡으로 사용합니다(상단 약 1% 꼬리).

    Returns
    -------
    WinsorizeExperimentResult

    Raises
    ------
    ValueError
        observations_by_variant가 비었거나, upper_tail이 허용값이 아니거나, 어느 실험그룹에도 유효 값이 없을 때.
    """
    #1) 인자 검증: 비어 있는 observations_by_variant·허용되지 않은 upper_tail을 조기에 거른다.
    if not observations_by_variant:
        raise ValueError("observations_by_variant must be non-empty.")
    if upper_tail is not None and upper_tail not in (0.05, 0.01):
        raise ValueError("upper_tail must be None, 0.05, or 0.01.")

    #2) 실험그룹별로 float 배열화 후 NaN 제거; 유효 관측이 없는 그룹은 풀을 만들 수 없으므로 예외.
    cleaned: Dict[str, np.ndarray] = {}
    for group_key, values in observations_by_variant.items():
        cleaned_arr = _strip_nan_float_array(values)
        if cleaned_arr.size == 0:
            raise ValueError(
                f"Group {group_key!r} has no valid (non-NaN) observations."
            )
        cleaned[group_key] = cleaned_arr

    #3) 모든 그룹의 유효 값을 이어 붙여 합동 풀(combined)을 만든다. 백분위·캡은 이 풀 기준으로만 계산한다.
    combined = np.concatenate(list(cleaned.values()))
    #4) 풀에서 90·95·99·100 백분위를 구하고, 로그/리포트용 문자열과 반환 필드에 쓴다.
    p90, p95, p99, p100 = np.percentile(combined, [90, 95, 99, 100])
    percentiles_str = (
        f"[90th: {p90:.2f}, 95th: {p95:.2f}, 99th: {p99:.2f}, 100th: {p100:.2f}]"
    )

    #5) upper_tail에 따라 상단 캡(upper_cap)만 정한다. None이면 이후 단계에서 복사만 한다.
    if upper_tail is None:
        upper_cap: Optional[float] = None
    elif upper_tail == 0.05:
        upper_cap = float(p95)
    else:
        upper_cap = float(p99)

    #6) 각 실험그룹 배열에 상한을 적용한다. 캡은 풀에서 한 번 정한 값을 모든 그룹에 동일하게 사용한다.
    winsorized: Dict[str, np.ndarray] = {}
    for group_key, cleaned_arr in cleaned.items():
        if upper_cap is None:
            winsorized[group_key] = cleaned_arr.copy()
        else:
            winsorized[group_key] = np.minimum(cleaned_arr, upper_cap)

    #7) 백분위·풀 크기·winsorize된 그룹별 배열을 묶어 반환한다.
    return WinsorizeExperimentResult(
        percentiles=percentiles_str,
        p90=float(p90),
        p95=float(p95),
        p99=float(p99),
        p100=float(p100),
        pool_n=int(combined.size),
        winsorized=winsorized,
    )


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
        한 행에 control_formula, treatment_formula, control_value, treatment_value,
        delta_relative, delta_absolute, p_value, CI_relative, CI_absolute,
        MSS_posthoc, statistic 컬럼.
    """
    epsilon = 1e-10

    if control_n <= 0 or treatment_n <= 0:
        raise ValueError("control_n and treatment_n must be positive.")
    if not (0 <= control_success <= control_n and 0 <= treatment_success <= treatment_n):
        raise ValueError("Success counts must be between 0 and respective n.")

    p1 = control_success / control_n
    p2 = treatment_success / treatment_n
    diff = p2 - p1

    # z-통계량·p-value: H0 하에서 합동추정량(p_hat)을 쓴 pooled 표준오차로 계산
    p_hat = (control_success + treatment_success) / (control_n + treatment_n)
    delta_se_pooled = np.sqrt(p_hat * (1 - p_hat) * (1 / control_n + 1 / treatment_n))
    if delta_se_pooled <= 0:
        delta_se_pooled = epsilon
    stat_z = diff / delta_se_pooled  # stat_z: z-통계량
    p_value = 2 * (1 - norm.cdf(abs(stat_z)))  # 양측검정 p-value

    # 신뢰구간: 비율 차이의 표준오차 (unpooled)
    delta_se_unpooled = np.sqrt(p1 * (1 - p1) / control_n + p2 * (1 - p2) / treatment_n)
    if delta_se_unpooled <= 0:
        delta_se_unpooled = epsilon

    # 대조군/실험군 지표 공식 및 값
    control_formula = f"{control_success}/{control_n}"
    treatment_formula = f"{treatment_success}/{treatment_n}"
    control_value = p1
    treatment_value = p2
    delta_absolute = diff
    delta_relative = (diff / p1) * 100 if p1 > epsilon else np.nan

    # 신뢰구간 (비율 차이: z 사용)
    z_value = norm.ppf(1 - alpha / 2)  # 양측검정 상위 α/2 지점
    delta_ci_lower = diff - z_value * delta_se_unpooled
    delta_ci_upper = diff + z_value * delta_se_unpooled
    CI_absolute = f"[{delta_ci_lower:.2f}, {delta_ci_upper:.2f}]"

    # 증감률(%) 신뢰구간 — 델타 메소드; 비율로 계산 후 CI_relative 표시 시 100 곱함
    if p1 > epsilon:
        uplift = diff / p1
        uplift_se = np.sqrt(
            (1 / p1**2) * p2 * (1 - p2) / treatment_n
            + (p2**2 / p1**4) * p1 * (1 - p1) / control_n
        )
        uplift_ci_lower = uplift - z_value * uplift_se
        uplift_ci_upper = uplift + z_value * uplift_se
        CI_relative = f"[{uplift_ci_lower * 100:.2f}%, {uplift_ci_upper * 100:.2f}%]"
    else:
        CI_relative = "[nan%, nan%]"

    # MSS_posthoc: Minimum Sample Size (최소 샘플 수) 대비 현재 비율 — 사후 계산
    beta = 1 - power  # beta: 제2종 오류 확률
    z_alpha = norm.ppf(1 - alpha / 2)  # z_alpha: 유의수준 α에 대한 z 임계값
    z_beta = norm.ppf(1 - beta)  # z_beta: 검정력(1-β)에 대한 z 값
    k = control_n / treatment_n if treatment_n > 0 else 0  # k: 대조군/실험군 샘플 비율

    if k <= 0 or not (0 < p1 < 1 and 0 < p2 < 1):
        MSS_posthoc = "0.0% (∞)"
    else:
        delta_abs = abs(p2 - p1)  # delta_abs: 비율 차이의 절댓값
        if delta_abs <= 0:
            MSS_posthoc = "0.0% (∞)"
        else:
            # n2_min: 검정력을 만족하기 위한 실험군 최소 샘플 수
            n2_min = ((z_alpha + z_beta) ** 2) * (
                (p1 * (1 - p1)) / k + p2 * (1 - p2)
            ) / (delta_abs**2)
            if np.isnan(n2_min) or n2_min <= 0:
                MSS_posthoc = "0.0% (∞)"
            else:
                min_n2 = int(np.ceil(n2_min))  # min_n2: 최소 샘플 수 (정수로 올림)
                ratio_pct = (treatment_n / min_n2) * 100  # ratio_pct: 현재 샘플 수 / 최소 샘플 수 (%)
                MSS_posthoc = f"{ratio_pct:.1f}% ({min_n2:,})"

    # 출력: control_value, treatment_value는 100 곱해 %로, delta_relative/ delta_absolute 포맷
    control_value_out = f"{control_value * 100:.2f}%"
    treatment_value_out = f"{treatment_value * 100:.2f}%"
    delta_relative_out = f"{delta_relative:.2f}%" if np.isfinite(delta_relative) else "nan%"
    delta_absolute_out = f"{delta_absolute * 100:.2f}%p"  # 비율 차이 → 퍼센트 포인트

    return pd.DataFrame(
        [
            {
                "control_formula": control_formula,
                "treatment_formula": treatment_formula,
                "control_value": control_value_out,
                "treatment_value": treatment_value_out,
                "delta_relative": delta_relative_out,
                "delta_absolute": delta_absolute_out,
                "p_value": round(p_value, 5),
                "CI_relative": CI_relative,
                "CI_absolute": CI_absolute,
                "MSS_posthoc": MSS_posthoc,
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
        한 행에 control_formula, treatment_formula, control_value, treatment_value,
        delta_relative, delta_absolute, p_value, CI_relative, CI_absolute,
        MSS_posthoc, statistic, df 컬럼.
    """
    x, y = _to_valid_arrays(control_values, treatment_values)
    n1, n2 = len(x), len(y)

    if n1 <= 1 or n2 <= 1:
        raise ValueError(
            "Each group must have at least 2 observations (for variance). "
            f"Got control n={n1}, treatment n={n2}."
        )

    # 리스트 기준: 관측수 = len, 누적값 = sum
    control_sum = float(np.sum(x))
    treatment_sum = float(np.sum(y))
    mu1 = float(np.mean(x))
    mu2 = float(np.mean(y))
    # 표본 분산 (ddof=1)
    s1_sq = float(np.var(x, ddof=1)) if n1 > 1 else 0.0
    s2_sq = float(np.var(y, ddof=1)) if n2 > 1 else 0.0

    # 대조군/실험군 지표 공식 (합계/관측수, 합은 소수 둘째 자리)
    control_formula = f"{control_sum:.2f}/{n1}"
    treatment_formula = f"{treatment_sum:.2f}/{n2}"
    control_value = mu1
    treatment_value = mu2

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
        t_value = t.ppf(1 - alpha / 2, df_welch)  # 양측검정 상위 α/2 지점
    else:
        p_value = np.nan
        t_value = np.nan

    diff = mu2 - mu1
    delta_absolute = diff
    epsilon = 1e-10
    delta_relative = (diff / mu1) * 100 if abs(mu1) > epsilon else np.nan

    # 신뢰구간 (평균 차이: t 사용)
    delta_ci_lower = diff - t_value * se
    delta_ci_upper = diff + t_value * se
    CI_absolute = f"[{delta_ci_lower:.2f}, {delta_ci_upper:.2f}]"

    # 증감률(%) 신뢰구간 — 델타 메소드; 비율로 계산 후 CI_relative 표시 시 100 곱함
    if abs(mu1) > epsilon:
        uplift = diff / mu1
        uplift_se = np.sqrt(
            (1 / mu1**2) * s2_sq / n2
            + (mu2**2 / mu1**4) * s1_sq / n1
        )
        uplift_ci_lower = uplift - t_value * uplift_se
        uplift_ci_upper = uplift + t_value * uplift_se
        CI_relative = f"[{uplift_ci_lower * 100:.2f}%, {uplift_ci_upper * 100:.2f}%]"
    else:
        CI_relative = "[nan%, nan%]"

    # MSS_posthoc: Minimum Sample Size (최소 샘플 수) - 평균 차이 기준 (Welch 가정, k = n1/n2) — 사후 계산
    beta = 1 - power  # beta: 제2종 오류 확률
    z_alpha = norm.ppf(1 - alpha / 2)  # z_alpha: 유의수준 α에 대한 z 임계값
    z_beta = norm.ppf(1 - beta)  # z_beta: 검정력(1-β)에 대한 z 값
    k = n1 / n2 if n2 > 0 else 0  # k: 대조군/실험군 샘플 비율

    if k <= 0:
        MSS_posthoc = "0.0% (∞)"
    else:
        delta_abs = abs(mu2 - mu1)  # delta_abs: 평균 차이의 절댓값
        if delta_abs <= 0 or np.isnan(delta_abs):
            MSS_posthoc = "0.0% (∞)"
        else:
            # n2_min: 검정력을 만족하기 위한 실험군 최소 샘플 수
            n2_min = ((z_alpha + z_beta) ** 2) * (s1_sq / k + s2_sq) / (delta_abs**2)
            if np.isnan(n2_min) or n2_min <= 0:
                MSS_posthoc = "0.0% (∞)"
            else:
                min_n2 = int(np.ceil(n2_min))  # min_n2: 최소 샘플 수 (정수로 올림)
                ratio_pct = (n2 / min_n2) * 100  # ratio_pct: 현재 샘플 수 / 최소 샘플 수 (%)
                MSS_posthoc = f"{ratio_pct:.1f}% ({min_n2:,})"

    # 출력: delta_relative 소수 둘째 자리 + '%', delta_absolute 소수 둘째 자리
    delta_relative_out = f"{delta_relative:.2f}%" if np.isfinite(delta_relative) else "nan%"
    delta_absolute_out = round(delta_absolute, 3)

    return pd.DataFrame(
        [
            {
                "control_formula": control_formula,
                "treatment_formula": treatment_formula,
                "control_value": control_value,
                "treatment_value": treatment_value,
                "delta_relative": delta_relative_out,
                "delta_absolute": delta_absolute_out,
                "p_value": round(p_value, 5),
                "CI_relative": CI_relative,
                "CI_absolute": CI_absolute,
                "MSS_posthoc": MSS_posthoc,
                "statistic": round(float(t_stat), 2),
                "df": round(float(df_welch), 2),
            }
        ]
    )
