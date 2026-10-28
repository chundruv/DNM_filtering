"""Data abstractions for genomic variant datasets."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, field
import hashlib

from .config import CHROM_COLUMN_CANDIDATES, POS_COLUMN_CANDIDATES


@dataclass
class VariantDataset:
    """
    Container for genomic variant data with metadata.
    
    This class provides a clean abstraction over raw DataFrames,
    handling variant-specific operations and caching.
    """
    variants: pd.DataFrame
    metadata: Dict[str, Any] = field(default_factory=dict)
    _hash: Optional[str] = None
    
    def __post_init__(self):
        """Optimize data types on initialization."""
        self._optimize_dtypes()
        
    def _optimize_dtypes(self):
        """Convert to optimal dtypes for memory and performance."""
        # Categorical columns for low cardinality
        categorical_cols = ['CHROM', 'var_type', 'REF', 'ALT']
        for col in categorical_cols:
            if col in self.variants.columns:
                self.variants[col] = self.variants[col].astype('category')
        
        # Downcast numeric types
        int_cols = self.variants.select_dtypes(include=['int64']).columns
        for col in int_cols:
            self.variants[col] = pd.to_numeric(self.variants[col], downcast='integer')
        
        float_cols = self.variants.select_dtypes(include=['float64']).columns
        for col in float_cols:
            self.variants[col] = pd.to_numeric(self.variants[col], downcast='float')
    
    @classmethod
    def from_path(
        cls,
        path: Path,
        **kwargs,
    ) -> 'VariantDataset':
        """Load from a file, dispatching on extension (.parquet vs TSV)."""
        path = Path(path)
        suffixes = {s.lower() for s in path.suffixes}
        if '.parquet' in suffixes or path.suffix.lower() == '.pq':
            return cls.from_parquet(path, **kwargs)
        # The TSV loader doesn't (yet) project columns, so drop the kwarg.
        kwargs.pop('columns', None)
        return cls.from_tsv(path, **kwargs)

    @classmethod
    def from_tsv(
        cls,
        path: Path,
        sample_col: str = 'SAMPLE',
        paternal_age_col: str = 'paternal_age',
        maternal_age_col: str = 'maternal_age',
        reference_col: str = 'ref',
        alternate_col: str = 'alt',
        required_cols: Optional[List[str]] = None,
        chrom_col: Optional[str] = None,
        pos_col: Optional[str] = None,
    ) -> 'VariantDataset':
        """Load from TSV file with automatic preprocessing."""
        df = pd.read_csv(path, sep='\t', on_bad_lines='skip', low_memory=False)

        # Clean quotes
        df.columns = df.columns.str.replace('"', '')
        for col in df.select_dtypes(include=['object']).columns:
            if df[col].astype(str).str.contains('"').any():
                df[col] = df[col].str.replace('"', '', regex=False)

        df = cls._normalize_columns(
            df, sample_col, paternal_age_col, maternal_age_col, reference_col, alternate_col,
            chrom_col, pos_col,
        )
        df = cls._add_computed_columns(df)

        # Dynamically detect and convert numeric columns
        # Try to convert all object-type columns to numeric
        object_cols = df.select_dtypes(include=['object']).columns
        for col in object_cols:
            # Count NaN values before conversion
            nan_before = df[col].isna().sum()

            # Try numeric conversion
            converted = pd.to_numeric(df[col], errors='coerce')
            nan_after = converted.isna().sum()

            # If conversion didn't create too many new NaNs, accept it
            # Allow up to 5% of values to be non-numeric
            threshold = len(df) * 0.05
            if nan_after - nan_before <= threshold:
                df[col] = converted

        # Ensure required_cols are converted if specified
        if required_cols:
            for col in required_cols:
                if col in df.columns and df[col].dtype == 'object':
                    df[col] = pd.to_numeric(df[col], errors='coerce')

        metadata = {
            'source': str(path),
            'n_variants': len(df),
            'variant_types': df['var_type'].value_counts().to_dict() if 'var_type' in df else {}
        }

        return cls(variants=df, metadata=metadata)

    @classmethod
    def from_parquet(
        cls,
        path: Path,
        sample_col: str = 'SAMPLE',
        paternal_age_col: str = 'paternal_age',
        maternal_age_col: str = 'maternal_age',
        reference_col: str = 'ref',
        alternate_col: str = 'alt',
        required_cols: Optional[List[str]] = None,
        chrom_col: Optional[str] = None,
        pos_col: Optional[str] = None,
        columns: Optional[List[str]] = None,
        dict_encode_strings: bool = True,
        downcast_numerics: bool = True,
    ) -> 'VariantDataset':
        """Load from a Parquet file.

        Much faster than from_tsv on large inputs because Parquet preserves dtypes,
        so we skip the per-column object->numeric scan that dominates TSV load time.

        Args:
            columns: If given, only load these columns. Names missing from the file
                are skipped (with a warning). Large memory saving on wide files.
            dict_encode_strings: Dictionary-encode string columns in Arrow before
                converting to pandas. Turns 200M-row SAMPLE strings into a Categorical
                with ~20K categories — usually a 10-50x peak-memory reduction.
            downcast_numerics: Cast float64->float32 unconditionally (precision loss
                accepted) and int64->int32 when all values fit. Halves numeric memory.
        """
        import pyarrow as pa
        import pyarrow.parquet as pq

        pf = pq.ParquetFile(str(path))
        available = set(pf.schema_arrow.names)

        if columns is not None:
            load_cols = [c for c in columns if c in available]
            missing = [c for c in columns if c not in available]
            if missing:
                import logging
                logging.getLogger(__name__).warning(
                    f"from_parquet: requested columns not in {path.name}, skipping: {missing}"
                )
        else:
            load_cols = None

        table = pq.read_table(str(path), columns=load_cols, use_threads=True)

        # Compact in Arrow before pandas materialization, since pandas object-dtype
        # strings + int64/float64 defaults are what blow up peak memory.
        if dict_encode_strings or downcast_numerics:
            new_arrays = []
            new_fields = []
            for i, field in enumerate(table.schema):
                col = table.column(i)
                t = field.type

                if dict_encode_strings and (
                    pa.types.is_string(t) or pa.types.is_large_string(t)
                ):
                    col = col.dictionary_encode()
                elif downcast_numerics and pa.types.is_floating(t) and t.bit_width > 32:
                    try:
                        col = col.cast(pa.float32(), safe=False)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        pass
                elif downcast_numerics and pa.types.is_integer(t) and t.bit_width > 32:
                    try:
                        col = col.cast(pa.int32(), safe=True)
                    except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
                        pass  # value out of int32 range, keep original

                new_arrays.append(col)
                new_fields.append(pa.field(field.name, col.type))

            table = pa.table(new_arrays, schema=pa.schema(new_fields))

        df = table.to_pandas(self_destruct=True)
        del table

        df = cls._normalize_columns(
            df, sample_col, paternal_age_col, maternal_age_col, reference_col, alternate_col,
            chrom_col, pos_col,
        )
        df = cls._add_computed_columns(df)

        if required_cols:
            for col in required_cols:
                if col in df.columns and df[col].dtype == 'object':
                    df[col] = pd.to_numeric(df[col], errors='coerce')

        metadata = {
            'source': str(path),
            'n_variants': len(df),
            'variant_types': df['var_type'].value_counts().to_dict() if 'var_type' in df else {}
        }

        return cls(variants=df, metadata=metadata)

    @staticmethod
    def _normalize_columns(
        df: pd.DataFrame,
        sample_col: str,
        paternal_age_col: str,
        maternal_age_col: str,
        reference_col: str,
        alternate_col: str,
        chrom_col: Optional[str] = None,
        pos_col: Optional[str] = None,
    ) -> pd.DataFrame:
        """Rename input columns to the canonical SAMPLE / paternal_age / maternal_age / REF / ALT / CHROM / POS names."""
        if sample_col != 'SAMPLE' and sample_col in df.columns:
            df.rename(columns={sample_col: 'SAMPLE'}, inplace=True)

        if paternal_age_col != 'paternal_age' and paternal_age_col in df.columns:
            df.rename(columns={paternal_age_col: 'paternal_age'}, inplace=True)

        if maternal_age_col != 'maternal_age' and maternal_age_col in df.columns:
            df.rename(columns={maternal_age_col: 'maternal_age'}, inplace=True)

        if reference_col != 'REF' and reference_col in df.columns:
            df.rename(columns={reference_col: 'REF'}, inplace=True)
        elif reference_col.lower() != 'ref':
            for old in ['Ref', 'Reference', 'ref', 'reference']:
                if old in df.columns and old != reference_col:
                    df.rename(columns={old: 'REF'}, inplace=True)
                    break

        if alternate_col != 'ALT' and alternate_col in df.columns:
            df.rename(columns={alternate_col: 'ALT'}, inplace=True)
        elif alternate_col.lower() != 'alt':
            for old in ['Alt', 'Alternate', 'alt', 'alternate', 'Variant']:
                if old in df.columns and old != alternate_col:
                    df.rename(columns={old: 'ALT'}, inplace=True)
                    break

        for canonical, configured, candidates in (
            ('CHROM', chrom_col, CHROM_COLUMN_CANDIDATES),
            ('POS', pos_col, POS_COLUMN_CANDIDATES),
        ):
            if configured and configured != canonical and configured in df.columns:
                df.drop(columns=[canonical], errors='ignore', inplace=True)
                df.rename(columns={configured: canonical}, inplace=True)
            elif canonical not in df.columns:
                for old in candidates:
                    if old in df.columns:
                        df.rename(columns={old: canonical}, inplace=True)
                        break

        return df

    @classmethod
    def _add_computed_columns(cls, df: pd.DataFrame) -> pd.DataFrame:
        """Add var_type, length, abs_length, midparage if the source columns are present."""
        if 'REF' in df.columns and 'ALT' in df.columns:
            df['var_type'] = cls._get_var_type(df)
            df['length'] = df['ALT'].str.len() - df['REF'].str.len()
            df['abs_length'] = np.abs(df['length'])

        if 'paternal_age' in df.columns and 'maternal_age' in df.columns:
            df['midparage'] = (df['paternal_age'] + df['maternal_age']) / 2

        return df
    
    @staticmethod
    def _get_var_type(df: pd.DataFrame) -> pd.Series:
        """Vectorized variant type detection."""
        ref_len = df['REF'].str.len()
        alt_len = df['ALT'].str.len()

        conditions = [
            (ref_len == 1) & (alt_len == 1),
            ref_len > alt_len,
            ref_len < alt_len
        ]
        choices = ['SNV', 'Deletion', 'Insertion']

        return pd.Series(np.select(conditions, choices, default='Other'), index=df.index)

    @staticmethod
    def _get_homopolymer_length(ref: str, alt: str) -> int:
        """
        Calculate the length of the longest single-base run in the longer allele.

        Scans the whole allele (REF for deletions, ALT for insertions), not just
        the inserted/deleted bases. Returns 1 for SNVs, 0 for non-string input.
        Not currently called by the package.

        Parameters
        ----------
        ref : str
            Reference allele sequence
        alt : str
            Alternate allele sequence

        Returns
        -------
        int
            Length of longest homopolymer run
        """
        if not isinstance(ref, str) or not isinstance(alt, str):
            return 0

        ref = ref.upper()
        alt = alt.upper()

        # For SNVs, homopolymer length is 1
        if len(ref) == 1 and len(alt) == 1:
            return 1

        # For indels, check the longer sequence (contains the repetitive region)
        seq_to_check = ref if len(ref) > len(alt) else alt

        if len(seq_to_check) == 0:
            return 0

        # Find longest run of repeated nucleotides
        max_length = 1
        current_length = 1
        current_base = seq_to_check[0]

        for base in seq_to_check[1:]:
            if base == current_base:
                current_length += 1
                max_length = max(max_length, current_length)
            else:
                current_base = base
                current_length = 1

        return max_length
    
    def filter_by_type(self, var_type: str) -> 'VariantDataset':
        """Return new dataset filtered by variant type."""
        filtered = self.variants[self.variants['var_type'] == var_type].copy()
        meta = self.metadata.copy()
        meta['filtered_type'] = var_type
        return VariantDataset(variants=filtered, metadata=meta)
    
    def apply_filters(self, params: Dict[str, Any]) -> 'VariantDataset':
        """Apply filtering parameters and return new dataset."""
        filtered = self.variants.copy()
        
        for param, value in params.items():
            if param.startswith('min_'):
                col = param[4:]  # Remove 'min_' prefix
                if col in filtered.columns:
                    filtered = filtered[filtered[col] >= value]
            elif param.startswith('max_'):
                col = param[4:]  # Remove 'max_' prefix  
                if col in filtered.columns:
                    filtered = filtered[filtered[col] <= value]
        
        return VariantDataset(variants=filtered, metadata=self.metadata)
    
    def count_by_sample(self, sample_col: str = 'SAMPLE') -> pd.Series:
        """Count variants per sample."""
        return self.variants.groupby(sample_col).size()
    
    def get_hash(self) -> str:
        """Get hash of dataset for caching."""
        if self._hash is None:
            # Hash based on shape and sample of data
            shape_str = f"{self.variants.shape}"
            sample_str = str(self.variants.sample(min(100, len(self.variants))).values)
            self._hash = hashlib.md5((shape_str + sample_str).encode()).hexdigest()
        return self._hash
    
    @property
    def summary(self) -> Dict[str, Any]:
        """Summary statistics."""
        return {
            'n_variants': len(self.variants),
            'n_samples': self.variants['SAMPLE'].nunique() if 'SAMPLE' in self.variants.columns else 0,
            'variant_types': self.variants['var_type'].value_counts().to_dict() if 'var_type' in self.variants.columns else {},
            'chromosomes': self.variants['CHROM'].nunique() if 'CHROM' in self.variants.columns else 0
        }
    
    def __len__(self):
        return len(self.variants)
    
    def __repr__(self):
        return f"VariantDataset(n_variants={len(self.variants)}, columns={list(self.variants.columns)})"
