from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


ATLAS_COLUMNS = ["chr", "start", "end", "startCpG", "endCpG", "target", "name", "direction"]
HG19_CPG_COUNT = 28_217_448


def load_marker_indices(atlas_path: Path):
    atlas = pd.read_csv(atlas_path, sep="\t", low_memory=False)
    atlas.columns = [str(name).strip().lstrip("#") for name in atlas.columns]

    if atlas.columns[:8].tolist() != ATLAS_COLUMNS:
        raise ValueError("Unexpected marker atlas columns")

    cell_types = list(atlas.columns[8:])
    if not cell_types:
        raise ValueError("No reference cell types in marker atlas")

    if atlas.duplicated(ATLAS_COLUMNS[:5] + ["target"]).any():
        raise ValueError("Duplicate marker-target entries")
    if not atlas["direction"].astype(str).str.upper().eq("U").all():
        raise ValueError("Marker direction must be U")

    for column in ("startCpG", "endCpG"):
        atlas[column] = pd.to_numeric(atlas[column], errors="raise").astype("int64")

    indices = {}
    counts = {}
    for cell_type in cell_types:
        markers = atlas.loc[atlas["target"].eq(cell_type), ["startCpG", "endCpG"]]
        if markers.empty:
            raise ValueError(f"Missing marker regions for {cell_type}")
        blocks = []
        for start, end in markers.itertuples(index=False, name=None):
            first = int(start) - 1
            last = int(end) - 1
            if first < 0 or last <= first:
                raise ValueError(f"Invalid CpG range for {cell_type}: {start}-{end}")
            blocks.append(np.arange(first, last, dtype=np.int64))
        indices[cell_type] = np.concatenate(blocks)
        counts[cell_type] = len(markers)

    return cell_types, indices, counts


def calculate_cfcpgratio(lbeta_path: Path, cell_types, marker_indices, marker_counts):
    size = lbeta_path.stat().st_size
    if not size or size % 4:
        raise ValueError("Invalid .lbeta file size")

    n_sites = size // 4
    lbeta = np.memmap(lbeta_path, dtype="<u2", mode="r", shape=(n_sites, 2))
    records = []
    try:
        for cell_type in cell_types:
            selected = marker_indices[cell_type]
            if int(selected.max()) >= n_sites:
                raise ValueError(f"CpG marker indices exceed .lbeta range for {cell_type}")

            methylated = np.asarray(lbeta[selected, 0], dtype=np.uint64)
            coverage = np.asarray(lbeta[selected, 1], dtype=np.uint64)
            if (methylated > coverage).any() or (coverage == 65535).any():
                raise ValueError(f"Invalid or saturated CpG counts for {cell_type}")

            total = int(coverage.sum(dtype=np.uint64))
            m = int(methylated.sum(dtype=np.uint64))
            u = total - m
            ratio = u / total if total else np.nan

            records.append({
                "CellType": cell_type,
                "Num_markers": marker_counts[cell_type],
                "CpG_slots": len(selected),
                "Covered_CpG_slots": int(np.count_nonzero(coverage)),
                "Total_CpG_observations": total,
                "Unmethylated_CpG_observations": u,
                "Methylated_CpG_observations": m,
                "cfCpGRatio_raw": ratio,
            })
    finally:
        del lbeta

    result = pd.DataFrame.from_records(records)
    values = result["cfCpGRatio_raw"].to_numpy(dtype=float)
    valid = np.isfinite(values)
    denominator = float(np.sum(values[valid]))
    result["cfCpGRatio_scaled"] = np.nan
    if denominator > 0:
        result["cfCpGRatio_scaled"] = 0.0
        result.loc[valid, "cfCpGRatio_scaled"] = values[valid] / denominator

    return result


def convert_pat(pat: Path, directory: Path, wgbstools: str, genome: str, threads: int):
    if not pat.name.endswith(".pat.gz"):
        raise ValueError("The PAT input must have the .pat.gz suffix")
    if not shutil.which(wgbstools):
        raise FileNotFoundError(f"Cannot find wgbstools executable: {wgbstools}")

    command = [
        wgbstools, "pat2beta", str(pat.resolve()),
        "--genome", genome, "--lbeta", "--threads", str(threads),
        "-o", str(directory),
    ]
    subprocess.run(command, check=True)

    stem = pat.name[:-len(".pat.gz")]
    lbeta = directory / f"{stem}.lbeta"
    if not lbeta.is_file():
        raise FileNotFoundError(f"wgbstools did not produce the expected file: {lbeta}")
    if genome == "hg19" and lbeta.stat().st_size != HG19_CPG_COUNT * 4:
        raise ValueError("The .lbeta file does not match the expected hg19 CpG count")
    return lbeta


def parse_args(argv=None):
    parser = argparse.ArgumentParser(prog="cfCpGRatio")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--pat", type=Path, help="Input .pat.gz file")
    inputs.add_argument("--lbeta", type=Path, help="Existing .lbeta file")
    parser.add_argument("--atlas", required=True, type=Path, help="Marker atlas TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output TSV")
    parser.add_argument("--genome", default="hg19", help="wgbstools genome name (default: hg19)")
    parser.add_argument("--wgbstools", default="wgbstools", help="Path to wgbstools executable")
    parser.add_argument("--threads", type=int, default=4, help="Threads for PAT conversion")
    parser.add_argument("--keep-lbeta", action="store_true", help="Keep converted .lbeta next to the output")
    return parser.parse_args(argv)


def run(args):
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if not args.atlas.is_file():
        raise FileNotFoundError(f"Marker atlas not found: {args.atlas}")
    input_path = args.pat if args.pat is not None else args.lbeta
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    cell_types, marker_indices, marker_counts = load_marker_indices(args.atlas)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.lbeta is not None:
        result = calculate_cfcpgratio(args.lbeta, cell_types, marker_indices, marker_counts)
    else:
        with tempfile.TemporaryDirectory(prefix="cfcpgratio-") as scratch:
            lbeta = convert_pat(args.pat, Path(scratch), args.wgbstools, args.genome, args.threads)
            result = calculate_cfcpgratio(lbeta, cell_types, marker_indices, marker_counts)
            if args.keep_lbeta:
                shutil.copy2(lbeta, args.output.parent / f"{args.output.stem}.lbeta")

    result.to_csv(args.output, sep="\t", index=False, na_rep="NA", float_format="%.12g")
    return result


def main(argv=None):
    try:
        args = parse_args(argv)
        result = run(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"cfCpGRatio: {error}", file=sys.stderr)
        return 1

    print(f"Processed {len(result)} cell types")
    print(f"Output: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
