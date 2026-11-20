from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ATLAS_KEY = ["chr", "start", "end", "startCpG", "endCpG"]
ATLAS_METADATA = ATLAS_KEY + ["target", "name", "direction"]
UINT16_BYTES_PER_CPG = 4


def load_marker_indices(atlas_path: Path) -> tuple[list[str], dict[str, np.ndarray], dict[str, int]]:
    atlas = pd.read_csv(atlas_path, sep="\t", low_memory=False)
    atlas.columns = [str(column).strip().lstrip("#") for column in atlas.columns]

    if atlas.columns[:8].tolist() != ATLAS_METADATA:
        raise ValueError(
            "Unexpected marker-atlas columns. Expected the first eight columns "
            f"to be: {', '.join(ATLAS_METADATA)}"
        )

    cell_types = list(atlas.columns[8:])
    if not cell_types:
        raise ValueError("No reference cell-type columns were found in the atlas.")

    required = set(ATLAS_METADATA)
    if not required.issubset(atlas.columns):
        raise ValueError(f"Marker atlas is missing columns: {sorted(required - set(atlas.columns))}")

    if not atlas["direction"].astype(str).str.upper().eq("U").all():
        raise ValueError("The supplied atlas contains markers that are not marked as U.")

    for column in ("startCpG", "endCpG"):
        atlas[column] = pd.to_numeric(atlas[column], errors="raise").astype("int64")

    marker_indices: dict[str, np.ndarray] = {}
    marker_counts: dict[str, int] = {}

    for cell_type in cell_types:
        rows = atlas.loc[atlas["target"].astype(str).eq(cell_type), ["startCpG", "endCpG"]]
        if rows.empty:
            raise ValueError(f"No marker regions were found for cell type: {cell_type}")

        blocks = []
        for start_cpg, end_cpg in rows.itertuples(index=False, name=None):
            start_index = int(start_cpg) - 1
            end_index = int(end_cpg) - 1

            if start_index < 0 or end_index <= start_index:
                raise ValueError(
                    f"Invalid CpG range for {cell_type}: {start_cpg}-{end_cpg}"
                )

            blocks.append(np.arange(start_index, end_index, dtype=np.int64))

        marker_indices[cell_type] = np.concatenate(blocks)
        marker_counts[cell_type] = len(rows)

    return cell_types, marker_indices, marker_counts


def calculate_cfcpgratio(
    lbeta_path: Path,
    cell_types: list[str],
    marker_indices: dict[str, np.ndarray],
    marker_counts: dict[str, int],
) -> pd.DataFrame:
    file_size = lbeta_path.stat().st_size
    if file_size == 0 or file_size % UINT16_BYTES_PER_CPG != 0:
        raise ValueError(
            "The .lbeta file size is invalid. Expected four bytes per CpG "
            "(two little-endian uint16 values: methylated and total coverage)."
        )

    n_cpg = file_size // UINT16_BYTES_PER_CPG
    lbeta = np.memmap(
        lbeta_path,
        mode="r",
        dtype="<u2",
        shape=(n_cpg, 2),
    )

    records = []
    try:
        for cell_type in cell_types:
            indices = marker_indices[cell_type]
            if int(indices.max()) >= n_cpg:
                raise ValueError(
                    f"Marker CpG index exceeds the .lbeta file for {cell_type}. "
                    "Check that the marker atlas and methylation data use the same "
                    "reference genome and CpG indexing."
                )

            methylated = np.asarray(lbeta[indices, 0], dtype=np.uint64)
            total = np.asarray(lbeta[indices, 1], dtype=np.uint64)

            if np.any(methylated > total):
                raise ValueError(
                    f"Methylated counts exceed total coverage for {cell_type}."
                )

            total_observations = int(total.sum(dtype=np.uint64))
            methylated_observations = int(methylated.sum(dtype=np.uint64))
            unmethylated_observations = total_observations - methylated_observations

            raw_ratio = (
                unmethylated_observations / total_observations
                if total_observations > 0
                else np.nan
            )

            records.append(
                {
                    "CellType": cell_type,
                    "Num_markers": marker_counts[cell_type],
                    "CpG_slots": int(len(indices)),
                    "Covered_CpG_slots": int(np.count_nonzero(total)),
                    "Total_CpG_observations": total_observations,
                    "Unmethylated_CpG_observations": unmethylated_observations,
                    "Methylated_CpG_observations": methylated_observations,
                    "cfCpGRatio_raw": raw_ratio,
                }
            )
    finally:
        del lbeta

    result = pd.DataFrame.from_records(records)
    finite = np.isfinite(result["cfCpGRatio_raw"].to_numpy(dtype=float))
    scale_total = float(result.loc[finite, "cfCpGRatio_raw"].sum())

    result["cfCpGRatio_scaled"] = np.nan
    if scale_total > 0:
        result["cfCpGRatio_scaled"] = 0.0
        result.loc[finite, "cfCpGRatio_scaled"] = (
            result.loc[finite, "cfCpGRatio_raw"] / scale_total
        )

    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate cfCpGRatio from a wgbstools .lbeta file."
    )
    parser.add_argument("--lbeta", required=True, type=Path, help="Input wgbstools .lbeta file.")
    parser.add_argument("--atlas", required=True, type=Path, help="Cell-type marker atlas in TSV format.")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV file.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    for path, label in ((args.lbeta, "lbeta file"), (args.atlas, "marker atlas")):
        if not path.is_file():
            print(f"Error: {label} not found: {path}", file=sys.stderr)
            return 2

    try:
        cell_types, marker_indices, marker_counts = load_marker_indices(args.atlas)
        result = calculate_cfcpgratio(
            args.lbeta,
            cell_types,
            marker_indices,
            marker_counts,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, sep="\t", index=False, na_rep="NA", float_format="%.12g")
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    print(f"Processed {len(cell_types)} reference cell types.")
    print(f"Results written to: {args.output}")
    print(
        "Sum of scaled estimates across cell types with available CpG coverage: "
        f"{result['cfCpGRatio_scaled'].sum(skipna=True):.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
