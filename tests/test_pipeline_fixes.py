import numpy as np
import optuna
import pandas as pd

from dnm_harmoniser.api import StageRunner
from dnm_harmoniser.config import (
    ColumnConfig,
    OptimisationConfig,
    PipelineConfig,
    Stage2Config,
    SymmetricRangeConstraint,
)
from dnm_harmoniser.data import VariantDataset
from dnm_harmoniser.pipeline import OptimisationPipeline


def _make_dataset(seed: int = 0, n_samples: int = 40) -> VariantDataset:
    """Synthetic SNV/insertion/deletion calls whose count rises with paternal age."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_samples):
        age = 20 + i % 25
        for _ in range(rng.poisson(5 + age)):
            ref, alt = rng.choice([("A", "T"), ("A", "AT"), ("AT", "A")], p=[0.8, 0.1, 0.1])
            rows.append({
                "sample_id": f"s{i}",
                "paternal_age": age,
                "maternal_age": age - 2,
                "ref": ref,
                "alt": alt,
                "quality": float(rng.uniform(0, 100)),
                "vaf": float(rng.uniform(0, 1)),
            })
    df = pd.DataFrame(rows)
    df = VariantDataset._normalize_columns(df, "sample_id", "paternal_age", "maternal_age", "ref", "alt")
    df = VariantDataset._add_computed_columns(df)
    return VariantDataset(variants=df)


def _config(**stage2) -> PipelineConfig:
    return PipelineConfig(
        stage2=Stage2Config(**stage2),
        optimisation=OptimisationConfig(
            columns=[
                ColumnConfig(name="quality", dtype="float", optimisation="maximum"),
                ColumnConfig(
                    name="vaf", dtype="float", optimisation="range",
                    range_constraint=SymmetricRangeConstraint(lower=0.1, scale=1.0),
                ),
            ]
        ),
        use_cache=False,
        max_workers=1,
    )


def test_outlier_removal_counts_samples_with_zero_filtered_dnms():
    data = _make_dataset()
    pipeline = OptimisationPipeline(_config(min_dnm_count=1, max_dnm_count=10_000))
    # Only s0's variants pass: every other sample ends up with 0 filtered DNMs
    s0_quality = data.variants.loc[data.variants["SAMPLE"] == "s0", "quality"].min()
    data.variants.loc[data.variants["SAMPLE"] != "s0", "quality"] = s0_quality - 1
    params = {t: {"min_quality": s0_quality} for t in ["SNV", "Insertion", "Deletion"]}

    keep, n_removed = pipeline._identify_outliers(data, params)
    keep_iter, n_removed_iter = pipeline._identify_outliers_with_params(data, params)

    n_samples = data.variants["SAMPLE"].nunique()
    assert keep == keep_iter == {"s0"}
    assert n_removed == n_removed_iter == n_samples - 1


def test_symmetric_range_never_suggests_min_above_max():
    data = _make_dataset()
    pipeline = OptimisationPipeline(_config())
    arrays_cols = {"quality": data.variants["quality"].to_numpy(), "vaf": data.variants["vaf"].to_numpy()}

    class _Arrays:
        filter_arrays = arrays_cols

    specs, linked = pipeline._build_search_space_specs(_Arrays(), "SNV")
    study = optuna.create_study(sampler=optuna.samplers.RandomSampler(seed=1))
    for _ in range(200):
        trial = study.ask()
        params = pipeline._suggest_params_for_arrays(trial, specs, linked, is_cmaes=False)
        assert 0.1 <= params["min_vaf"] <= 0.5
        assert params["min_vaf"] <= params["max_vaf"]
        study.tell(trial, 0.0)


def test_stage_runner_runs_all_stages():
    data = _make_dataset(seed=1)
    reference = _make_dataset(seed=2)
    runner = StageRunner(_config(min_dnm_count=1, max_dnm_count=10_000))

    targets = runner.calculate_targets(reference)
    warmup = runner.warmup(data, targets, n_trials=30)
    clean = runner.remove_outliers(data, warmup)
    final = runner.optimize(clean, targets, n_trials=20)

    assert set(final) == {"SNV", "Insertion", "Deletion"}
    assert all("min_quality" in p for p in final.values())
    assert set(clean.variants["SAMPLE"]) <= set(data.variants["SAMPLE"])
    # Per-call threshold overrides don't leak into the config
    runner.remove_outliers(data, warmup, min_dnm=3, max_dnm=50)
    assert (runner.config.stage2.min_dnm_count, runner.config.stage2.max_dnm_count) == (1, 10_000)
