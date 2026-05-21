"""Stream a (possibly gzipped) TSV into a Parquet file.

Reads with pyarrow's multi-threaded CSV reader and writes batches incrementally,
so peak memory stays near a single block instead of loading the full file.
"""

from pathlib import Path
from typing import List, Optional
import logging
import time

import typer

logger = logging.getLogger(__name__)


def convert_tsv_to_parquet(
    input_path: Path,
    output_path: Path,
    columns: Optional[List[str]] = None,
    block_size_mb: int = 64,
    compression: str = "snappy",
    skip_bad_rows: bool = True,
) -> int:
    """Stream-convert a TSV to Parquet. Returns the number of rows written."""
    import pyarrow.csv as pacsv
    import pyarrow.parquet as pq

    read_opts = pacsv.ReadOptions(
        block_size=block_size_mb * 1024 * 1024,
        use_threads=True,
    )
    parse_opts = pacsv.ParseOptions(
        delimiter="\t",
        invalid_row_handler=(lambda row: "skip") if skip_bad_rows else None,
    )
    convert_opts = pacsv.ConvertOptions(
        include_columns=columns,
        strings_can_be_null=True,
        null_values=["", "NA", "NaN", "nan", "NULL", "null", "."],
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    reader = pacsv.open_csv(
        str(input_path),
        read_options=read_opts,
        parse_options=parse_opts,
        convert_options=convert_opts,
    )

    writer: Optional[pq.ParquetWriter] = None
    n_rows = 0
    n_batches = 0
    t0 = time.time()

    try:
        while True:
            try:
                batch = reader.read_next_batch()
            except StopIteration:
                break

            if writer is None:
                writer = pq.ParquetWriter(
                    str(output_path),
                    batch.schema,
                    compression=compression,
                )
            writer.write_batch(batch)
            n_rows += batch.num_rows
            n_batches += 1

            if n_batches % 20 == 0:
                elapsed = time.time() - t0
                rate = n_rows / elapsed if elapsed > 0 else 0
                logger.info(f"  {n_batches} batches, {n_rows:,} rows ({rate:,.0f} rows/s)")
    finally:
        if writer is not None:
            writer.close()
        reader.close()

    elapsed = time.time() - t0
    logger.info(f"Wrote {n_rows:,} rows in {n_batches} batches in {elapsed:.1f}s")
    return n_rows


def convert_command(
    input: Path = typer.Argument(..., help="Input TSV (.tsv or .tsv.gz)"),
    output: Path = typer.Argument(..., help="Output Parquet file"),
    columns: Optional[List[str]] = typer.Option(
        None, "--columns", "-c",
        help="Restrict to these columns (repeatable). Default: keep all."
    ),
    block_size: int = typer.Option(64, "--block-size-mb", help="CSV read block size, MB"),
    compression: str = typer.Option(
        "snappy", "--compression",
        help="Parquet compression: snappy (fast), zstd (smaller), none"
    ),
    keep_bad_rows: bool = typer.Option(
        False, "--keep-bad-rows",
        help="Fail on malformed rows instead of skipping (default: skip, matching the TSV loader)"
    ),
):
    """Convert a TSV to Parquet for fast reloads.

    The output keeps the original column names, so existing config files work unchanged
    once the loader points at the .parquet file.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )

    if not input.exists():
        typer.echo(f"Error: input file not found: {input}", err=True)
        raise typer.Exit(1)

    typer.echo(f"Converting {input} -> {output}")
    if columns:
        typer.echo(f"  Keeping columns: {', '.join(columns)}")

    n_rows = convert_tsv_to_parquet(
        input_path=input,
        output_path=output,
        columns=columns,
        block_size_mb=block_size,
        compression=compression,
        skip_bad_rows=not keep_bad_rows,
    )

    size_mb = output.stat().st_size / 1024 / 1024
    typer.echo(f"Done: {n_rows:,} rows, {size_mb:,.1f} MB at {output}")
