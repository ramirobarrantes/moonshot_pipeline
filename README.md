# Somatic Variant Calling + Copy Number Analysis Pipeline

A Nextflow DSL2 pipeline for **somatic variant calling** (MuSE), **VEP annotation**, and **copy number alteration estimation** (ichorCNA) from matched tumor/normal CRAM pairs.

## Overview

```
CRAM (tumor + normal)
    │
    └── samtools/convert ──► BAM
            │
            ├── MuSE call ──► MuSE sump ──► Somatic VCF
            │                                     │
            │                               VEP annotation ──► Annotated VCF
            │
            └── hmmcopy/readcounter ──► ichorCNA ──► CNA segments + plots
                (tumor + normal)
```

| Step | Tool | Purpose |
|------|------|---------|
| 1 | samtools convert | CRAM → BAM for tumor and normal |
| 2a | MuSE call | Pre-filter somatic variant positions |
| 2b | MuSE sump | Apply tier cutoffs → somatic VCF |
| 2c | Ensembl VEP | Annotate somatic VCF (offline cache) |
| 3a | hmmcopy readcounter | Count reads in genomic windows (WIG) |
| 3b | ichorCNA run | Estimate tumor fraction and copy number segments |

All modules follow [nf-core/modules](https://github.com/nf-core/modules) conventions.

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
    --vep_cache /path/to/.vep \
    --outdir results \
    -resume
```

For **WXS (exome)** data, add `--wgs false` or use `-profile wxs`.

VEP annotation is skipped if `--vep_cache` is omitted.

---

## Samplesheet Format

The pipeline expects a CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `patient` | Unique sample identifier (used as the output prefix) |
| `tumor_cram` | Absolute path to the tumor CRAM file |
| `tumor_crai` | Absolute path to the tumor CRAM index (.crai) |
| `normal_cram` | Absolute path to the matched normal CRAM file |
| `normal_crai` | Absolute path to the matched normal CRAM index (.crai) |

**Example** (`assets/samplesheet.csv`):

```csv
patient,tumor_cram,tumor_crai,normal_cram,normal_crai
MS008_TUMOR,/data/MS008_TUMOR.recal.cram,/data/MS008_TUMOR.recal.cram.crai,/data/MS008_NORMAL.recal.cram,/data/MS008_NORMAL.recal.cram.crai
MS008_TUMOR2,/data/MS008_TUMOR2.recal.cram,/data/MS008_TUMOR2.recal.cram.crai,/data/MS008_NORMAL.recal.cram,/data/MS008_NORMAL.recal.cram.crai
```

Multiple tumors can share the same normal — the pipeline handles the pairing correctly by joining on `patient` id.

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
| `--gc_wig` | GC content WIG for ichorCNA (must match bin size) |
| `--map_wig` | Mappability WIG for ichorCNA (must match bin size) |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wgs` | `true` | WGS mode (`-G`). Set `false` for WXS (`-E`). |
| `--vep_cache` | — | Path to VEP offline cache directory (`~/.vep`) |
| `--vep_genome` | `GRCh38` | Assembly name passed to VEP `--assembly` |
| `--vep_cache_version` | `114` | VEP cache version |
| `--centromere` | — | Centromere positions text file for ichorCNA |
| `--panel_of_normals` | — | Panel of normals RDS for ichorCNA |
| `--rep_time_wig` | — | Replication timing WIG for ichorCNA |
| `--exons` | — | Exon BED for ichorCNA annotation |
| `--outdir` | `./results` | Output directory |

---

## Reference Files

### ichorCNA WIG files (`--gc_wig`, `--map_wig`)

WIG files containing per-window GC content and mappability scores that correct systematic read-depth biases before copy number estimation. **The bin size must match the readcounter window** (pipeline default: 1 Mb).

Pre-built files for hg38 at 1 Mb resolution:
- `gc_hg38_1000kb.wig`
- `map_hg38_1000kb.wig`

Available from the ichorCNA repository:
https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata

For other bin sizes, generate custom files with `gcCounter` and `mapCounter` from the HMMcopy suite.

### dbSNP (`--dbsnp`, `--dbsnp_tbi`)

Required by MuSE sump. Must be bgzipped and tabix-indexed. Download from NCBI:
https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/

### VEP cache (`--vep_cache`)

Download the GRCh38 cache with:

```bash
vep_install -a cf -s homo_sapiens -y GRCh38 --CACHE_VERSION 114 -c /path/to/.vep
```

Or use the Docker/Apptainer image to run the installer:

```bash
apptainer exec docker://ensemblorg/ensembl-vep:release_114.0 \
    vep_install -a cf -s homo_sapiens -y GRCh38 --CACHE_VERSION 114 -c ~/.vep
```

---

## Output Structure

```
results/
├── bam/
│   ├── tumor/           # Converted tumor BAMs (.bam + .bai)
│   └── normal/          # Converted normal BAMs (.bam + .bai)
├── muse/
│   ├── call/            # MuSE intermediate files (*.MuSE.txt)
│   └── sump/            # Somatic VCFs (*.vcf.gz + .tbi)
├── vep/                 # Annotated VCFs (*.vep.vcf.gz + .tbi + .summary.html)
├── hmmcopy/
│   ├── tumor/           # Tumor read count WIG files
│   └── normal/          # Normal read count WIG files
└── ichorcna/
    ├── *.seg            # Copy number segments
    ├── *.cna.seg        # CNA segment file
    ├── *.seg.txt        # Segment text table
    ├── *.params.txt     # Estimated tumor fraction and ploidy
    ├── *.correctedDepth.txt
    ├── *.RData          # R objects for reanalysis
    └── plots/           # Genome-wide CNA plots (PDF)
```

---

## Resource Configuration

Default resource labels (`nextflow.config`):

| Label | CPUs | Memory | Time |
|-------|------|--------|------|
| `process_single` | 1 | 1 GB | 30 min |
| `process_low` | 2 | 6 GB | 2 h |
| `process_medium` | 6 | 24 GB | 6 h |
| `process_high` | 12 | 64 GB | 12 h |

Global limits: 72 CPUs, 512 GB RAM, 168 h per job.

Override per-process or globally with a custom config:

```bash
nextflow run main.nf ... -c my_cluster.config
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| (default) | Slurm + Apptainer, WGS mode |
| `local` | Local executor (no Slurm); containers still pulled via Wave |
| `wxs` | WXS/exome mode — MuSE sump uses `-E` instead of `-G` |

---

## Requirements

- **Nextflow** ≥ 24.04.0
- **Apptainer** (Singularity) — containers pulled automatically via Wave
- **Slurm** — default executor; use `-profile local` to run without it

---

## License

MIT
