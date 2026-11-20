# cfCpGRatio

cfCpGRatio estimates  cell-type-specific cfDNA signals directly from CpG methylation observations in selected marker regions. 
## Installation

Python 3.9 or newer is required. A Conda environment is recommended:

```bash
conda create -n cfCpGRatio python=3.11 pip -y
conda activate cfCpGRatio
conda install -c conda-forge -c bioconda numpy pandas scipy samtools htslib -y
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

The external dependencies required by `wgbstools` (including samtools, bgzip and tabix) are listed in its installation instructions. If wgbstools is already installed on the analysis server, use the existing executable instead.

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


