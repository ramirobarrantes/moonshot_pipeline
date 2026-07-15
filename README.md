# Somatic Variant Calling + Copy Number Analysis Pipeline

A Nextflow DSL2 pipeline for **somatic variant calling** (MuSE) and **copy number alteration estimation** (ichorCNA) from matched tumor/normal CRAM pairs.

## Overview

```
CRAM (tumor + normal)
    │
    ├── samtools/convert ──► BAM
    │       │
    │       ├── MuSE call ──► MuSE sump ──► Somatic VCF
    │       │
    │       └── hmmcopy/readcounter ──► ichorCNA ──► CNA segments + plots
    │
    └── (matched normal also feeds ichorCNA as baseline)
```

| Step | Tool | Purpose |
|------|------|---------|
| 1 | samtools convert | CRAM → BAM conversion for tumor and normal |
| 2a | MuSE call | Pre-filter somatic variant positions |
| 2b | MuSE sump | Apply tier cutoffs → final somatic VCF |
| 3a | hmmcopy readcounter | Count reads in genomic windows (WIG) |
| 3b | ichorCNA run | Estimate tumor fraction and copy number segments |

All modules are sourced from [nf-core/modules](https://github.com/nf-core/modules).

---

## Quick Start

```bash
nextflow run main.nf \
    --input assets/samplesheet.csv \
    --fasta /path/to/genome.fa \
    --fai /path/to/genome.fa.fai \
    --dbsnp /path/to/dbsnp.vcf.gz \
    --dbsnp_tbi /path/to/dbsnp.vcf.gz.tbi \
    --gc_wig /path/to/gc_hg38_1000kb.wig \
    --map_wig /path/to/map_hg38_1000kb.wig \
    --outdir results \
    -resume
```

For **WXS (exome)** data, add `--wgs false` or use `-profile wxs`.

---

## Samplesheet Format

The pipeline expects a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `patient` | Unique patient/sample identifier |
| `tumor_cram` | Path to the tumor CRAM file |
| `tumor_crai` | Path to the tumor CRAM index (.crai) |
| `normal_cram` | Path to the matched normal CRAM file |
| `normal_crai` | Path to the matched normal CRAM index (.crai) |

**Example** (`assets/samplesheet.csv`):

```csv
patient,tumor_cram,tumor_crai,normal_cram,normal_crai
PATIENT_01,/data/tumor_01.cram,/data/tumor_01.cram.crai,/data/normal_01.cram,/data/normal_01.cram.crai
PATIENT_02,/data/tumor_02.cram,/data/tumor_02.cram.crai,/data/normal_02.cram,/data/normal_02.cram.crai
```

Each row represents one tumor/normal pair. The `patient` value is used as the sample ID throughout the pipeline.

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to the samplesheet CSV |
| `--fasta` | Reference genome FASTA |
| `--fai` | Reference genome FASTA index (.fai) |
| `--dbsnp` | dbSNP VCF file (bgzipped) |
| `--dbsnp_tbi` | dbSNP VCF tabix index (.tbi) |
| `--gc_wig` | GC content WIG file for ichorCNA (matching bin size) |
| `--map_wig` | Mappability WIG file for ichorCNA (matching bin size) |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wgs` | `true` | WGS mode (`-G` for MuSE sump). Set to `false` for WXS (`-E`). |
| `--centromere` | — | Centromere positions text file for ichorCNA |
| `--panel_of_normals` | — | Panel of normals RDS file for ichorCNA |
| `--rep_time_wig` | — | Replication timing WIG file for ichorCNA |
| `--exons` | — | Exon BED file for ichorCNA annotation |
| `--outdir` | `./results` | Output directory |

---

## Reference Files

### ichorCNA reference files (`--gc_wig` and `--map_wig`)

These are WIG-format files containing per-window GC content and mappability scores. They correct systematic biases in read depth before copy number estimation.

- **`gc_wig`** — GC content per genomic bin. Corrects for GC bias in library preparation.
- **`map_wig`** — Mappability per genomic bin. Filters unreliable regions (repeats, centromeres).

**The bin size must match the readcounter window size** (default: 1,000,000 bp / 1 Mb).

Pre-built files for common genomes are available from the ichorCNA repository:
https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata

For hg38 at 1 Mb resolution:
- `gc_hg38_1000kb.wig`
- `map_hg38_1000kb.wig`

For other bin sizes (500kb, 100kb, 50kb, 10kb), use the corresponding files or generate custom ones with `gcCounter` and `mapCounter` from the HMMcopy suite.

### dbSNP (`--dbsnp` and `--dbsnp_tbi`)

MuSE sump requires a dbSNP VCF to annotate variant positions. Download from NCBI:
https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/

Ensure the VCF is bgzipped and indexed with tabix.

---

## Output Structure

```
results/
├── bam/
│   ├── tumor/          # Converted tumor BAMs
│   └── normal/         # Converted normal BAMs
├── muse/
│   ├── call/           # MuSE intermediate call files
│   └── sump/           # Final somatic VCFs (.vcf.gz + .tbi)
├── hmmcopy/
│   ├── tumor/          # Tumor read count WIG files
│   └── normal/         # Normal read count WIG files
└── ichorcna/
    ├── *.seg           # Copy number segments
    ├── *.params.txt    # Estimated tumor fraction and ploidy
    ├── *.RData         # R objects for reanalysis
    └── plots/          # CNA genome-wide plots (PDF)
```

---

## Resource Configuration

Default resource labels are defined in `nextflow.config`:

| Label | CPUs | Memory | Time |
|-------|------|--------|------|
| `process_single` | 1 | 1 GB | 30 min |
| `process_low` | 2 | 6 GB | 2 h |
| `process_medium` | 6 | 24 GB | 6 h |
| `process_high` | 12 | 64 GB | 12 h |

Override in a custom config or on the command line:

```bash
nextflow run main.nf --input samplesheet.csv ... \
    -c my_cluster.config
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| (default) | WGS mode, Docker + Wave enabled |
| `wxs` | WXS/exome mode (MuSE sump uses `-E` instead of `-G`) |

---

## Requirements

- **Nextflow** ≥ 24.04.0
- **Docker** (or Singularity/Podman with appropriate config changes)
- **Wave** is enabled by default for container resolution

---

## License

MIT
