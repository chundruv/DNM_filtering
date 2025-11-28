"""Data abstractions for genomic variant datasets."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, List, Any
from dataclasses import dataclass, field
import hashlib


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
    def from_tsv(
        cls,
        path: Path,
        sample_col: str = 'SAMPLE',
        paternal_age_col: str = 'paternal_age',
        maternal_age_col: str = 'maternal_age',
        reference_col: str = 'ref',
        alternate_col: str = 'alt',
        required_cols: Optional[List[str]] = None,
        max_length: int = 20
    ) -> 'VariantDataset':
        """Load from TSV file with automatic preprocessing."""
        df = pd.read_csv(path, sep='\t', on_bad_lines='skip', low_memory=False)

        # Clean quotes
        df.columns = df.columns.str.replace('"', '')
        for col in df.select_dtypes(include=['object']).columns:
            if df[col].astype(str).str.contains('"').any():
                df[col] = df[col].str.replace('"', '', regex=False)

        # Standardize column names
        if sample_col != 'SAMPLE' and sample_col in df.columns:
            df.rename(columns={sample_col: 'SAMPLE'}, inplace=True)

        if paternal_age_col != 'paternal_age' and paternal_age_col in df.columns:
            df.rename(columns={paternal_age_col: 'paternal_age'}, inplace=True)

        if maternal_age_col != 'maternal_age' and maternal_age_col in df.columns:
            df.rename(columns={maternal_age_col: 'maternal_age'}, inplace=True)

        # Standardize reference and alternate column names
        if reference_col != 'REF' and reference_col in df.columns:
            df.rename(columns={reference_col: 'REF'}, inplace=True)
        elif reference_col.lower() != 'ref':
            # Handle common variations if not already standardized
            for old in ['Ref', 'Reference', 'ref', 'reference']:
                if old in df.columns and old != reference_col:
                    df.rename(columns={old: 'REF'}, inplace=True)
                    break

        if alternate_col != 'ALT' and alternate_col in df.columns:
            df.rename(columns={alternate_col: 'ALT'}, inplace=True)
        elif alternate_col.lower() != 'alt':
            # Handle common variations if not already standardized
            for old in ['Alt', 'Alternate', 'alt', 'alternate', 'Variant']:
                if old in df.columns and old != alternate_col:
                    df.rename(columns={old: 'ALT'}, inplace=True)
                    break
        
        # Calculate variant type and derived metrics
        if 'REF' in df.columns and 'ALT' in df.columns:
            df['var_type'] = cls._get_var_type(df)
            df['length'] = df['ALT'].str.len() - df['REF'].str.len()
            df['abs_length'] = np.abs(df['length'])  # Absolute indel length for filtering
            df['homopolymer_length'] = df.apply(lambda row: cls._get_homopolymer_length(row['REF'], row['ALT']), axis=1)
            df = df[df['abs_length'] < max_length]
        
        # Calculate midparage if parental ages exist
        if 'paternal_age' in df.columns and 'maternal_age' in df.columns:
            df['midparage'] = (df['paternal_age'] + df['maternal_age']) / 2
        
        # Convert numeric columns (including common variant filtering columns)
        numeric_cols = [
            'paternal_age', 'maternal_age', 'VAF', 'QUAL', 'DP',
            'IMF', 'DNM', 'MQ', 'DeNovoCNN_prob', 'nparAADn0',
            'child_coverage', 'father_coverage', 'mother_coverage',
            'FS', 'GQ', 'PL', 'AD', 'IDV', 'RPBZ', 'MQBZ', 'BQBZ', 'MQSBZ', 'SCBZ', 'SGB', 'MQ0F'
        ]
        if required_cols:
            numeric_cols.extend([c for c in required_cols if c not in numeric_cols])

        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        metadata = {
            'source': str(path),
            'n_variants': len(df),
            'variant_types': df['var_type'].value_counts().to_dict() if 'var_type' in df else {}
        }
        
        return cls(variants=df, metadata=metadata)
    
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
        Calculate the length of the longest homopolymer (repeated nucleotide) in REF or ALT.

        For deletions: checks the deleted sequence (REF - ALT)
        For insertions: checks the inserted sequence (ALT - REF)
        For SNVs: returns 1 (single nucleotide)

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
