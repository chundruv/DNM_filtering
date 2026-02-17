"""VAF distribution quality metrics for DNM filtering.

Provides four metrics to quantify VAF distribution quality, targeting
the left-skew that indicates false positive de novo calls. No reference
VAF distribution is needed — metrics use theoretical expectations for
true heterozygous germline de novo variants.

Metrics:
    1. Mean VAF — pulled down by heavy left tail
    2. Low-VAF proportion — fraction below a threshold
    3. Skewness — third standardised moment, captures asymmetry
    4. KL divergence — against a truncated normal, captures full shape

All metrics are variant-type-specific, accounting for the natural
left skew of indels due to alignment bias.
"""

import numpy as np
from scipy import stats
from typing import Dict, Optional, Tuple
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Expected values for true heterozygous germline de novos.
# Derived empirically from clean datasets (NHS-GMS, AoU patterns).
# ------------------------------------------------------------------

EXPECTED_MEAN: Dict[str, float] = {
    "SNV": 0.48,
    "Insertion": 0.44,
    "Deletion": 0.46,
}

EXPECTED_DISTRIBUTION_PARAMS: Dict[str, Dict[str, float]] = {
    "SNV": {"mu": 0.50, "sigma": 0.08},
    "Insertion": {"mu": 0.46, "sigma": 0.09},
    "Deletion": {"mu": 0.48, "sigma": 0.09},
}

LOW_VAF_THRESHOLDS: Dict[str, float] = {
    "SNV": 0.35,
    "Insertion": 0.30,
    "Deletion": 0.32,
}

EXPECTED_SKEWNESS: Dict[str, float] = {
    "SNV": -0.1,
    "Insertion": -0.3,
    "Deletion": -0.2,
}


@dataclass
class VAFQualityResult:
    """Container for VAF quality metrics for a single variant type."""

    variant_type: str
    n_variants: int
    mean_vaf: float
    low_vaf_proportion: float
    skewness: float
    kl_divergence: float

    # Quality scores normalised to [0, 1] where 1 = best quality
    mean_score: float = 0.0
    proportion_score: float = 0.0
    skewness_score: float = 0.0
    kl_score: float = 0.0
    overall_score: float = 0.0

    low_vaf_threshold: float = 0.35
    expected_mean: float = 0.48

    def to_dict(self) -> Dict:
        return {
            "variant_type": self.variant_type,
            "n_variants": self.n_variants,
            "mean_vaf": round(self.mean_vaf, 4),
            "low_vaf_proportion": round(self.low_vaf_proportion, 4),
            "skewness": round(self.skewness, 4),
            "kl_divergence": round(self.kl_divergence, 4),
            "mean_score": round(self.mean_score, 3),
            "proportion_score": round(self.proportion_score, 3),
            "skewness_score": round(self.skewness_score, 3),
            "kl_score": round(self.kl_score, 3),
            "overall_score": round(self.overall_score, 3),
        }

    def summary(self) -> str:
        return (
            f"{self.variant_type} (n={self.n_variants:,}): "
            f"mean={self.mean_vaf:.3f} (expected {self.expected_mean:.2f}), "
            f"low_proportion={self.low_vaf_proportion:.3f} "
            f"(below {self.low_vaf_threshold:.2f}), "
            f"skewness={self.skewness:.3f}, "
            f"KL_div={self.kl_divergence:.4f}, "
            f"overall_quality={self.overall_score:.2f}"
        )


# ------------------------------------------------------------------
# Fast metric functions for use inside optimisation objective.
# These operate on raw numpy arrays for minimal overhead per trial.
# ------------------------------------------------------------------


def compute_vaf_penalty(
    vaf: np.ndarray,
    variant_type: str,
    metric: str = "kl_divergence",
    vaf_min: float = 0.25,
    vaf_max: float = 0.75,
    n_bins: int = 50,
) -> float:
    """Compute a single VAF quality penalty for use in the objective function.

    This is the main entry point for the pipeline's objective function.
    Returns a penalty value where 0 = perfect quality and higher = worse.

    Parameters
    ----------
    vaf : np.ndarray
        VAF values of filtered variants (already within [vaf_min, vaf_max]).
    variant_type : str
        One of 'SNV', 'Insertion', 'Deletion'.
    metric : str
        Which metric to use: 'mean', 'proportion', 'skewness', 'kl_divergence'.
    vaf_min, vaf_max : float
        VAF bounds used in filtering.
    n_bins : int
        Number of bins for KL divergence histogram.

    Returns
    -------
    float
        Penalty value >= 0. Lower = better quality.
    """
    if len(vaf) < 10:
        # Too few variants to assess quality meaningfully
        return 0.0

    if metric == "mean":
        return _penalty_mean(vaf, variant_type)
    elif metric == "proportion":
        return _penalty_proportion(vaf, variant_type)
    elif metric == "skewness":
        return _penalty_skewness(vaf, variant_type)
    elif metric == "kl_divergence":
        return _penalty_kl(vaf, variant_type, vaf_min, vaf_max, n_bins)
    else:
        raise ValueError(f"Unknown VAF quality metric: {metric}. "
                         f"Choose from: mean, proportion, skewness, kl_divergence")


def _penalty_mean(vaf: np.ndarray, variant_type: str) -> float:
    """Penalty based on deviation of mean VAF from expected.

    Normalised so that a 0.05 deviation from expected gives penalty ~1.0.
    """
    expected = EXPECTED_MEAN.get(variant_type, 0.48)
    deviation = abs(np.mean(vaf) - expected)
    # Scale: 0.05 deviation -> penalty ~1.0
    return (deviation / 0.05) ** 2


def _penalty_proportion(vaf: np.ndarray, variant_type: str) -> float:
    """Penalty based on excess low-VAF proportion beyond expected baseline.

    Only penalises the excess above what's normal for the variant type.
    Normalised so that 10% excess proportion gives penalty ~1.0.
    """
    threshold = LOW_VAF_THRESHOLDS.get(variant_type, 0.35)
    proportion = float(np.mean(vaf < threshold))

    # Expected baseline in clean data
    expected_baseline = {"SNV": 0.03, "Insertion": 0.07, "Deletion": 0.05}
    baseline = expected_baseline.get(variant_type, 0.05)

    excess = max(0, proportion - baseline)
    # Scale: 0.10 excess -> penalty ~1.0
    return (excess / 0.10) ** 2


def _penalty_skewness(vaf: np.ndarray, variant_type: str) -> float:
    """Penalty based on deviation of skewness from expected.

    Penalises deviation in either direction — bimodal contamination
    can produce positive skew despite a heavy left tail.
    Normalised so that 0.5 units deviation gives penalty ~1.0.
    """
    expected = EXPECTED_SKEWNESS.get(variant_type, -0.1)
    skewness = float(stats.skew(vaf, bias=False))
    deviation = abs(skewness - expected)
    # Scale: 0.5 deviation -> penalty ~1.0
    return (deviation / 0.5) ** 2


def _penalty_kl(
    vaf: np.ndarray,
    variant_type: str,
    vaf_min: float = 0.25,
    vaf_max: float = 0.75,
    n_bins: int = 50,
) -> float:
    """Penalty based on KL divergence from theoretical truncated normal.

    Returns the KL divergence directly — already a natural penalty
    where 0 = perfect match. Typical clean data: 0.01-0.02.
    Typical noisy data: 0.1-0.5+.

    Normalised so KL ~0.1 gives penalty ~1.0.
    """
    params = EXPECTED_DISTRIBUTION_PARAMS.get(
        variant_type, {"mu": 0.48, "sigma": 0.09}
    )
    mu = params["mu"]
    sigma = params["sigma"]

    bin_edges = np.linspace(vaf_min, vaf_max, n_bins + 1)

    obs_counts, _ = np.histogram(vaf, bins=bin_edges)
    obs_counts = obs_counts.astype(np.float64) + 1.0
    obs_probs = obs_counts / obs_counts.sum()

    a = (vaf_min - mu) / sigma
    b = (vaf_max - mu) / sigma
    trunc_dist = stats.truncnorm(a, b, loc=mu, scale=sigma)
    ref_probs = np.diff(trunc_dist.cdf(bin_edges))
    ref_probs = ref_probs + 1e-10
    ref_probs = ref_probs / ref_probs.sum()

    kl = float(stats.entropy(obs_probs, ref_probs))
    # Scale: KL 0.1 -> penalty ~1.0
    return (kl / 0.1) ** 2


# ------------------------------------------------------------------
# Full-featured class for standalone analysis / reporting
# ------------------------------------------------------------------


class VAFQualityMetrics:
    """Compute VAF distribution quality metrics.

    Parameters
    ----------
    vaf_min, vaf_max : float
        VAF bounds used in filtering.
    n_bins : int
        Number of bins for KL divergence.
    expected_means, low_vaf_thresholds, distribution_params, expected_skewness : dict
        Override defaults per variant type.
    """

    def __init__(
        self,
        vaf_min: float = 0.25,
        vaf_max: float = 0.75,
        n_bins: int = 50,
        expected_means: Optional[Dict[str, float]] = None,
        low_vaf_thresholds: Optional[Dict[str, float]] = None,
        distribution_params: Optional[Dict[str, Dict[str, float]]] = None,
        expected_skewness: Optional[Dict[str, float]] = None,
    ):
        self.vaf_min = vaf_min
        self.vaf_max = vaf_max
        self.n_bins = n_bins
        self.expected_means = {**EXPECTED_MEAN, **(expected_means or {})}
        self.low_vaf_thresholds = {**LOW_VAF_THRESHOLDS, **(low_vaf_thresholds or {})}
        self.distribution_params = {
            **EXPECTED_DISTRIBUTION_PARAMS, **(distribution_params or {})
        }
        self.expected_skewness = {**EXPECTED_SKEWNESS, **(expected_skewness or {})}

    def compute_all(
        self,
        vaf: np.ndarray,
        variant_type: str = "SNV",
        weights: Optional[Dict[str, float]] = None,
    ) -> VAFQualityResult:
        """Compute all metrics and quality scores for a VAF array."""
        vaf = np.asarray(vaf, dtype=np.float64)

        if len(vaf) == 0:
            logger.warning(f"No variants for {variant_type}")
            return VAFQualityResult(
                variant_type=variant_type, n_variants=0,
                mean_vaf=np.nan, low_vaf_proportion=np.nan,
                skewness=np.nan, kl_divergence=np.nan,
            )

        expected_mean = self.expected_means.get(variant_type, 0.48)
        threshold = self.low_vaf_thresholds.get(variant_type, 0.35)

        mean_vaf = float(np.mean(vaf))
        low_prop = float(np.mean(vaf < threshold))
        skewness = float(stats.skew(vaf, bias=False))
        kl_div = self._compute_kl(vaf, variant_type)

        # Quality scores [0, 1], 1 = best
        mean_score = float(np.exp(-0.5 * ((mean_vaf - expected_mean) / 0.03) ** 2))

        expected_baseline = {"SNV": 0.03, "Insertion": 0.07, "Deletion": 0.05}
        excess = max(0, low_prop - expected_baseline.get(variant_type, 0.05))
        prop_score = float(np.exp(-10.0 * excess))

        expected_skew = self.expected_skewness.get(variant_type, -0.1)
        skew_score = float(np.exp(-2.5 * abs(skewness - expected_skew)))

        kl_score = float(np.exp(-10.0 * kl_div))

        w = weights or {"mean": 1.0, "proportion": 1.0, "skewness": 1.0, "kl": 1.0}
        total_weight = sum(w.values())
        overall = (
            w.get("mean", 1.0) * mean_score
            + w.get("proportion", 1.0) * prop_score
            + w.get("skewness", 1.0) * skew_score
            + w.get("kl", 1.0) * kl_score
        ) / total_weight

        return VAFQualityResult(
            variant_type=variant_type, n_variants=len(vaf),
            mean_vaf=mean_vaf, low_vaf_proportion=low_prop,
            skewness=skewness, kl_divergence=kl_div,
            mean_score=mean_score, proportion_score=prop_score,
            skewness_score=skew_score, kl_score=kl_score,
            overall_score=overall,
            low_vaf_threshold=threshold, expected_mean=expected_mean,
        )

    def _compute_kl(self, vaf: np.ndarray, variant_type: str) -> float:
        params = self.distribution_params.get(variant_type, {"mu": 0.48, "sigma": 0.09})
        mu, sigma = params["mu"], params["sigma"]
        bin_edges = np.linspace(self.vaf_min, self.vaf_max, self.n_bins + 1)

        obs_counts, _ = np.histogram(vaf, bins=bin_edges)
        obs_counts = obs_counts.astype(np.float64) + 1.0
        obs_probs = obs_counts / obs_counts.sum()

        a = (self.vaf_min - mu) / sigma
        b = (self.vaf_max - mu) / sigma
        trunc_dist = stats.truncnorm(a, b, loc=mu, scale=sigma)
        ref_probs = np.diff(trunc_dist.cdf(bin_edges)) + 1e-10
        ref_probs = ref_probs / ref_probs.sum()

        return float(stats.entropy(obs_probs, ref_probs))

    def summarise_cohort(
        self, df, vaf_col: str = "proband_VAF", type_col: str = "variant_type",
        variant_types: Optional[list] = None,
    ) -> Dict[str, VAFQualityResult]:
        """Compute metrics for all variant types in a DataFrame."""
        if variant_types is None:
            variant_types = sorted(df[type_col].unique())
        results = {}
        for vtype in variant_types:
            vaf = df.loc[df[type_col] == vtype, vaf_col].dropna().values
            results[vtype] = self.compute_all(vaf, variant_type=vtype)
            logger.info(results[vtype].summary())
        return results
