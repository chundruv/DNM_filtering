from pathlib import Path

import pandas as pd

from dnm_harmoniser.config import PipelineConfig
from dnm_harmoniser.data import VariantDataset


def test_parquet_projection_preserves_coordinate_columns(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            'sample_id': ['s1', 's2'],
            'paternal_age': [40, 42],
            'maternal_age': [35, 36],
            'ref': ['A', 'G'],
            'alt': ['T', 'GA'],
            'chr': ['1', '2'],
            'position': [101, 202],
        }
    )
    parquet_path = tmp_path / 'variants.parquet'
    df.to_parquet(parquet_path, index=False)

    config = PipelineConfig()
    dataset = VariantDataset.from_path(
        parquet_path,
        columns=config.optimisation.get_input_columns('data'),
    )

    assert 'CHROM' in dataset.variants.columns
    assert 'POS' in dataset.variants.columns
    assert dataset.variants['CHROM'].astype(str).tolist() == ['1', '2']
    assert dataset.variants['POS'].tolist() == [101, 202]


def test_parquet_projection_keeps_configured_coordinate_columns(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            'sample_id': ['s1', 's2'],
            'paternal_age': [40, 42],
            'maternal_age': [35, 36],
            'ref': ['A', 'G'],
            'alt': ['T', 'GA'],
            'Chr': ['chr1', 'chr2'],
            'start': [101, 202],
            'unused': [0, 1],
        }
    )
    parquet_path = tmp_path / 'variants.parquet'
    df.to_parquet(parquet_path, index=False)

    config = PipelineConfig()
    config.optimisation.position_column_data = 'start'
    dataset = VariantDataset.from_path(
        parquet_path,
        pos_col=config.optimisation.position_column_data,
        columns=config.optimisation.get_input_columns('data'),
    )

    # 'Chr' is auto-detected; 'start' only via the config
    assert dataset.variants['CHROM'].astype(str).tolist() == ['chr1', 'chr2']
    assert dataset.variants['POS'].tolist() == [101, 202]
    assert 'unused' not in dataset.variants.columns
