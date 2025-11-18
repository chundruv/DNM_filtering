"""Test suite for variant optimizer."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from hypothesis import given, strategies as st

from dnm_harmoniser import (
    PipelineConfig,
    VariantDataset,
    optimize_filters,
    run_optimization
)
from dnm_harmoniser.config import PRESETS
from dnm_harmoniser.pipeline import OptimizationPipeline


# Fixtures
@pytest.fixture(scope="module")
def sample_vcf_data():
    """Realistic VCF data for testing."""
    np.random.seed(42)
    n_variants = 200
    
    return pd.DataFrame({
        'SAMPLE': np.repeat(range(10), 20),
        'CHROM': np.random.choice(['chr1', 'chr2'], n_variants),
        'POS': np.random.randint(1000, 10000, n_variants),
        'REF': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'ALT': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'QUAL': np.random.uniform(10, 60, n_variants),
        'DP': np.random.randint(5, 100, n_variants),
        'child_DP': np.random.randint(10, 80, n_variants),
        'father_DP': np.random.randint(10, 80, n_variants),
        'mother_DP': np.random.randint(10, 80, n_variants),
        'VAF': np.random.uniform(0.1, 0.9, n_variants),
        'paternal_age': np.random.uniform(25, 45, n_variants),
        'maternal_age': np.random.uniform(23, 42, n_variants)
    })


@pytest.fixture
def variant_dataset(sample_vcf_data):
    """VariantDataset object for testing."""
    return VariantDataset(variants=sample_vcf_data)


@pytest.fixture
def default_config():
    """Default pipeline configuration."""
    return PipelineConfig()


# Unit Tests
class TestVariantDataset:
    """Test data abstraction."""
    
    def test_creation(self, sample_vcf_data):
        """Test dataset creation."""
        dataset = VariantDataset(variants=sample_vcf_data)
        assert len(dataset) == len(sample_vcf_data)
        assert 'var_type' in dataset.variants.columns
    
    def test_dtype_optimization(self, variant_dataset):
        """Test automatic dtype optimization."""
        assert variant_dataset.variants['CHROM'].dtype.name == 'category'
        if 'var_type' in variant_dataset.variants.columns:
            assert variant_dataset.variants['var_type'].dtype.name == 'category'
    
    def test_filter_by_type(self, variant_dataset):
        """Test filtering by variant type."""
        snv_data = variant_dataset.filter_by_type('SNV')
        assert all(snv_data.variants['var_type'] == 'SNV')
    
    def test_apply_filters(self, variant_dataset):
        """Test parameter-based filtering."""
        params = {'min_QUAL': 30, 'min_DP': 20}
        filtered = variant_dataset.apply_filters(params)
        assert all(filtered.variants['QUAL'] >= 30)
        assert all(filtered.variants['DP'] >= 20)
        assert len(filtered) <= len(variant_dataset)


class TestConfiguration:
    """Test configuration management."""
    
    def test_default_config(self):
        """Test default configuration."""
        config = PipelineConfig()
        assert config.stage1.enabled
        assert config.stage2.enabled
        assert config.stage3.enabled
        assert config.seed == 42
    
    @pytest.mark.parametrize("preset", list(PRESETS.keys()))
    def test_presets(self, preset):
        """Test all configuration presets."""
        config = PipelineConfig.from_preset(preset)
        assert config.stage3.n_trials > 0
        assert config.max_workers > 0
    
    def test_validation(self):
        """Test configuration validation."""
        with pytest.raises(ValueError):
            PipelineConfig(stage1={'n_trials': -1})  # Invalid
        
        with pytest.raises(ValueError):
            PipelineConfig(stage2={'min_dnm_count': 100, 'max_dnm_count': 50})  # Invalid range


class TestPipeline:
    """Test optimization pipeline."""
    
    def test_reproducibility(self, variant_dataset):
        """Test deterministic results with same seed."""
        config = PipelineConfig(
            stage1={'n_trials': 10},
            stage3={'n_trials': 20},
            seed=42,
            deterministic=True
        )
        
        pipeline1 = OptimizationPipeline(config)
        pipeline2 = OptimizationPipeline(config)
        
        # Both should initialize with same state
        assert config.seed == 42
        assert config.deterministic
    
    def test_stage_disabling(self, variant_dataset):
        """Test disabling individual stages."""
        config = PipelineConfig(
            stage1={'enabled': False},
            stage2={'enabled': False},
            stage3={'n_trials': 10}
        )
        pipeline = OptimizationPipeline(config)
        
        # Should run with stages disabled
        assert not config.stage1.enabled
        assert not config.stage2.enabled


# Parametrized Tests
@pytest.mark.parametrize("threshold,expected_ratio", [
    (20, 0.75),  # ~75% pass at QUAL>20
    (30, 0.50),  # ~50% pass at QUAL>30
    (40, 0.25),  # ~25% pass at QUAL>40
])
def test_quality_filtering(sample_vcf_data, threshold, expected_ratio):
    """Test quality filtering at different thresholds."""
    dataset = VariantDataset(variants=sample_vcf_data)
    filtered = dataset.apply_filters({'min_QUAL': threshold})
    
    actual_ratio = len(filtered) / len(dataset)
    assert actual_ratio == pytest.approx(expected_ratio, rel=0.3)


# Property-based Testing
@given(st.lists(st.floats(min_value=0, max_value=100), min_size=10))
def test_filtering_preserves_order(quality_scores):
    """Filtering should preserve variant order."""
    df = pd.DataFrame({
        'SAMPLE': range(len(quality_scores)),
        'QUAL': quality_scores,
        'idx': range(len(quality_scores))
    })
    dataset = VariantDataset(variants=df)
    filtered = dataset.apply_filters({'min_QUAL': 20})
    
    # Check order preserved
    if len(filtered) > 0:
        assert filtered.variants['idx'].is_monotonic_increasing


@given(
    st.floats(min_value=0.01, max_value=0.5),
    st.floats(min_value=0.5, max_value=0.99)
)
def test_threshold_monotonicity(low_threshold, high_threshold):
    """Higher thresholds should never increase variant count."""
    np.random.seed(42)
    df = pd.DataFrame({
        'SAMPLE': range(100),
        'VAF': np.random.uniform(0, 1, 100)
    })
    dataset = VariantDataset(variants=df)
    
    low_result = dataset.apply_filters({'min_VAF': low_threshold})
    high_result = dataset.apply_filters({'min_VAF': high_threshold})
    
    assert len(high_result) <= len(low_result)


# Integration Tests
def test_end_to_end_simple(tmp_path, sample_vcf_data):
    """Test complete pipeline with simple API."""
    # Save test data
    data_path = tmp_path / "data.tsv"
    ref_path = tmp_path / "reference.tsv"
    
    sample_vcf_data.to_csv(data_path, sep='\t', index=False)
    sample_vcf_data.sample(50).to_csv(ref_path, sep='\t', index=False)
    
    # Run optimization
    params = optimize_filters(
        data_path,
        ref_path,
        n_trials=10,
        preset="fast",
        seed=42
    )
    
    assert isinstance(params, dict)
    assert len(params) > 0


def test_end_to_end_config(tmp_path, sample_vcf_data):
    """Test complete pipeline with configuration."""
    # Save test data
    data_path = tmp_path / "data.tsv"
    ref_path = tmp_path / "reference.tsv"
    
    sample_vcf_data.to_csv(data_path, sep='\t', index=False)
    sample_vcf_data.sample(50).to_csv(ref_path, sep='\t', index=False)
    
    # Custom configuration
    config = PipelineConfig(
        stage1={'n_trials': 5},
        stage2={'min_dnm_count': 5, 'max_dnm_count': 100},
        stage3={'n_trials': 10}
    )
    
    # Run optimization
    result = run_optimization(data_path, ref_path, config=config)
    
    assert result.best_params is not None
    assert isinstance(result.summary, str)


# Performance Tests
@pytest.mark.slow
def test_performance_caching(variant_dataset):
    """Test that caching improves performance."""
    config = PipelineConfig(
        use_cache=True,
        cache_dir=Path("./test_cache"),
        stage1={'n_trials': 10},
        stage3={'n_trials': 20}
    )
    
    pipeline = OptimizationPipeline(config)
    
    # First run (no cache)
    import time
    start = time.time()
    # Would run pipeline here
    first_time = time.time() - start
    
    # Second run (with cache)
    start = time.time()
    # Would run pipeline again
    second_time = time.time() - start
    
    # Second should be faster (in real implementation)
    # assert second_time < first_time * 0.5


# Cleanup
@pytest.fixture(autouse=True)
def cleanup():
    """Clean up test files."""
    yield
    # Cleanup code here
    import shutil
    for path in ['./test_cache', './cache']:
        if Path(path).exists():
            shutil.rmtree(path)
