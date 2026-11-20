# cfCpGRatio

cfCpGRatio estimates tissue- and cell-type-associated cfDNA signals from CpG methylation observations in selected marker regions. Its primary input is a `wgbstools` methylation-pattern file (`.pat.gz`).

For cell type *i*, the raw estimate is

\[ \mathrm{cfCpGRatio}_i = \frac{U_i}{U_i+M_i} \]

where *U* and *M* denote the total unmethylated and methylated CpG observations across the markers assigned to that cell type. The program also reports estimates normalized to sum to one across the reference cell types within each sample. These estimates are not calibrated absolute tissue fractions.

## Installation

Python 3.9 or newer is required. A Conda environment is recommended:

```bash
conda create -n cfCpGRatio python=3.11 pip -y
conda activate cfCpGRatio
conda install -c conda-forge -c bioconda numpy pandas scipy samtools htslib tabix -y
```

Install `wgbstools` according to its [official instructions](https://github.com/nloyfer/wgbs_tools), and configure the reference genome used with the atlas. For hg19:

```bash
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py
export PATH="$PWD:$PATH"
wgbstools init_genome hg19
cd ..
```

`wgbstools` requires `samtools`, `tabix`, and `bgzip`. Install these tools with Conda before running cfCpGRatio; they are not installed by `pip install .` or `pip install -r requirements.txt`. The `tabix` package is listed explicitly because installing `htslib` alone may not make the `tabix` executable available in every Conda environment. `bgzip` is provided by HTSlib.

Confirm that the executables are available in the active environment:

```bash
command -v samtools
command -v tabix
command -v bgzip
```

All three commands should print a path. If `tabix` is missing from an existing environment, install it without recreating the environment:

```bash
conda install -c conda-forge -c bioconda tabix -y
```

If `wgbstools` is already installed on the analysis server, use the existing executable instead of reinstalling it.

Clone this repository and install the command-line program:

```bash
git clone https://github.com/liqiang-hku/cfCpGRatio.git
cd cfCpGRatio
python -m pip install .
cfCpGRatio --help
```

Alternatively, install Python dependencies without installing the command-line entry point:

```bash
python -m pip install -r requirements.txt
python src/cfcpgratio.py --help
```

## Input

The main input is a bgzip-compressed methylation-pattern file produced with `wgbstools`, ending in `.pat.gz`. The conversion step requires an hg19 genome initialized in `wgbstools` and the corresponding marker atlas.

`wgbstools pat2beta --lbeta` converts the PAT input to CpG-level methylated and total counts. cfCpGRatio reads those counts and calculates the ratio for each cell type. The temporary `.lbeta` file is removed after processing unless `--keep-lbeta` is specified.

The marker atlas must contain `chr`, `start`, `end`, `startCpG`, `endCpG`, `target`, `name`, and `direction` as its first eight columns, followed by reference cell-type columns. The `target` column indicates the cell type assigned to each marker region. Marker direction must be `U`. Use the marker coordinates and CpG indices that match the input PAT genome.

## Usage

```bash
cfCpGRatio \
  --pat sample.pat.gz \
  --atlas markers/Atlas.U250.l4.hg19.full.tsv \
  --genome hg19 \
  --output results/sample_cfCpGRatio.tsv \
  --threads 4
```

If `wgbstools` is not in the shell's `PATH`, provide its executable path:

```bash
cfCpGRatio \
  --pat sample.pat.gz \
  --atlas markers/Atlas.U250.l4.hg19.full.tsv \
  --wgbstools /path/to/wgbstools \
  --output results/sample_cfCpGRatio.tsv
```

An existing `.lbeta` file may also be analyzed directly:

```bash
cfCpGRatio \
  --lbeta sample.lbeta \
  --atlas markers/Atlas.U250.l4.hg19.full.tsv \
  --output results/sample_cfCpGRatio.tsv
```

## Output

One tab-separated row is written per reference cell type. The output includes CpG coverage totals, unmethylated and methylated observation counts, `cfCpGRatio_raw`, and `cfCpGRatio_scaled`.

The raw ratio is `NA` when no CpG observations cover a cell type's marker regions. For scaled composition estimates, a missing raw ratio contributes zero if the normalization denominator is positive. When all ratios are undefined or sum to zero, the scaled estimates remain `NA`.

## Reference

Loyfer N, et al. A DNA methylation atlas of normal human cell types. *Nature*. 2023;613:355–364. https://doi.org/10.1038/s41586-022-05580-6

Loyfer N, et al. wgbstools: a computational suite for DNA methylation sequencing data analysis. *Life Science Alliance*. 2026. See https://github.com/nloyfer/wgbs_tools

## Code and data

This repository contains the cfCpGRatio calculation. It does not include raw clinical sequencing files or the external marker atlas. Public datasets and analysis inputs are described in the accompanying article.
