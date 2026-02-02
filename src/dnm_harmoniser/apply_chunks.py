import pandas as pd
import sys

def get_homopolymer_length(ref: str, alt: str) -> int:
    """..."""  # same as before
    if not isinstance(ref, str) or not isinstance(alt, str):
        return 0
    ref, alt = ref.upper(), alt.upper()
    if len(ref) == 1 and len(alt) == 1:
        return 1
    seq = ref if len(ref) > len(alt) else alt
    if len(seq) == 0:
        return 0
    max_len = current_len = 1
    current_base = seq[0]
    for base in seq[1:]:
        if base == current_base:
            current_len += 1
            max_len = max(max_len, current_len)
        else:
            current_base = base
            current_len = 1
    return max_len

def process_chunk(chunk, ref_col, alt_col):
    """Filter a chunk by length and homopolymer criteria."""
    chunk = chunk.copy()
    chunk['length'] = chunk[alt_col].str.len() - chunk[ref_col].str.len()
    chunk = chunk[chunk['length'].abs() < 20]
    
    if len(chunk) == 0:
        return chunk
    
    chunk['homopolymer_length'] = [
        get_homopolymer_length(r, a) 
        for r, a in zip(chunk[ref_col], chunk[alt_col])
    ]
    return chunk[chunk['homopolymer_length'] < 8]

def process_file_chunked(input_path, output_path, ref_col, alt_col, chunksize=50_000):
    """Process a file in chunks, writing incrementally."""
    header_written = False
    total_in = total_out = 0
    
    for chunk in pd.read_csv(input_path, sep='\t', chunksize=chunksize, dtype=str):
        total_in += len(chunk)
        filtered = process_chunk(chunk, ref_col, alt_col)
        
        if len(filtered) > 0:
            filtered.to_csv(
                output_path, 
                sep='\t', 
                index=False, 
                mode='a', 
                header=not header_written
            )
            header_written = True
            total_out += len(filtered)
    
    print(f"{input_path}: {total_in} → {total_out} variants")

# Process both files
process_file_chunked(
    sys.argv[1],
    sys.argv[2],
    ref_col=sys.argv[3], alt_col=sys.argv[4]
)