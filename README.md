# Somatic Variant Calling + Copy Number Analysis Pipeline

A Nextflow DSL2 pipeline for **somatic variant calling** (MuSE + Sentieon TNScope), **VEP annotation**, and **copy number alteration estimation** (ichorCNA) from matched tumor/normal CRAM pairs.

## Overview

```
CRAM (tumor + normal)
    │
    ├── samtools/convert ──► BAM
    │       │
    │       ├── MuSE call ──► MuSE sump ──────────────────┐
    │       │                                              │
    │       └── hmmcopy/readcounter ──► ichorCNA           ├──► VEP annotation ──► Annotated VCFs
    │           (tumor + normal)       CNA + tumor fraction│
    │                                                      │
    └── Sentieon TNhaplotyper2 ──► TNfilter ───────────────┘
        (CRAMs direct, no BAM needed)
```

| Step | Tool | Purpose |
|------|------|---------|
| 1 | samtools convert | CRAM → BAM for tumor and normal |
| 2a | MuSE call | Pre-filter somatic variant positions |
| 2b | MuSE sump | Apply tier cutoffs → somatic VCF |
| 3a | hmmcopy readcounter | Count reads in genomic windows (WIG) |
| 3b | ichorCNA run | Estimate tumor fraction and copy number segments |
| 4a | Sentieon TNhaplotyper2 | Somatic calling + orientation bias + contamination model |
| 4b | Sentieon TNfilter | Filter raw VCF → final somatic VCF |
| 5 | Ensembl VEP | Annotate MuSE and TNScope VCFs (offline cache) |

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
    --germline_vcf /path/to/af-only-gnomad.hg38.vcf.gz \
    --germline_vcf_tbi /path/to/af-only-gnomad.hg38.vcf.gz.tbi \
    --contamination_vcf /path/to/small_exac_common_3.hg38.vcf.gz \
    --contamination_vcf_tbi /path/to/small_exac_common_3.hg38.vcf.gz.tbi \
    --vep_cache /path/to/.vep \
    --outdir results \
    -resume
```

For **WXS (exome)** data, add `--wgs false` or use `-profile wxs`.

VEP is skipped if `--vep_cache` is omitted.

---

## Samplesheet Format

The pipeline expects a CSV with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `patient` | yes | Unique sample identifier (used as the output prefix) |
| `tumor_cram` | yes | Absolute path to the tumor CRAM file |
| `tumor_crai` | yes | Absolute path to the tumor CRAM index (.crai) |
| `normal_cram` | yes | Absolute path to the matched normal CRAM file |
| `normal_crai` | yes | Absolute path to the matched normal CRAM index (.crai) |
| `tumor_sample` | TNScope | SM tag in the tumor CRAM `@RG` header — must match exactly |
| `normal_sample` | TNScope | SM tag in the normal CRAM `@RG` header — must match exactly |

`tumor_sample` and `normal_sample` **must exactly match the `SM:` tags** in the CRAM `@RG` header lines — Sentieon will error if they don't. To look them up:

```bash
samtools view -H your.cram | grep "^@RG"
# e.g. @RG  ID:...  SM:MS008_TUMOR  ...
```

**Example** (`assets/samplesheet.csv`):

```csv
patient,tumor_cram,tumor_crai,normal_cram,normal_crai,tumor_sample,normal_sample
MS008_TUMOR,/data/MS008_TUMOR.recal.cram,/data/MS008_TUMOR.recal.cram.crai,/data/MS008_NORMAL.recal.cram,/data/MS008_NORMAL.recal.cram.crai,MS008_TUMOR,MS008_NORMAL
MS008_TUMOR2,/data/MS008_TUMOR2.recal.cram,/data/MS008_TUMOR2.recal.cram.crai,/data/MS008_NORMAL.recal.cram,/data/MS008_NORMAL.recal.cram.crai,MS008_TUMOR2,MS008_NORMAL
```

Multiple tumors can share the same normal — the pipeline pairs them correctly by joining on `patient` id.

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to the samplesheet CSV |
| `--fasta` | Reference genome FASTA |
| `--fai` | Reference genome FASTA index (.fai) |
| `--dbsnp` | dbSNP VCF (bgzipped, for MuSE sump) |
| `--dbsnp_tbi` | dbSNP VCF tabix index (.tbi) |
| `--gc_wig` | GC content WIG for ichorCNA (must match bin size) |
| `--map_wig` | Mappability WIG for ichorCNA (must match bin size) |

### TNScope (Sentieon) — required to run Steps 4–5

| Parameter | Description |
|-----------|-------------|
| `--germline_vcf` | Germline population VCF (e.g. `af-only-gnomad.hg38.vcf.gz`) |
| `--germline_vcf_tbi` | Tabix index for germline VCF |
| `--contamination_vcf` | Common variant VCF for contamination model (e.g. `small_exac_common_3.hg38.vcf.gz`) |
| `--contamination_vcf_tbi` | Tabix index for contamination VCF |
### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--wgs` | `true` | WGS mode (`-G`). Set `false` for WXS (`-E`). |
| `--vep_cache` | — | Path to VEP offline cache directory |
| `--vep_genome` | `GRCh38` | Assembly name passed to VEP `--assembly` |
| `--vep_cache_version` | `114` | VEP cache version |
| `--centromere` | — | Centromere positions file for ichorCNA |
| `--panel_of_normals` | — | Panel of normals RDS for ichorCNA |
| `--rep_time_wig` | — | Replication timing WIG for ichorCNA |
| `--exons` | — | Exon BED for ichorCNA annotation |
| `--outdir` | `./results` | Output directory |

---

## Reference Files

### ichorCNA WIG files (`--gc_wig`, `--map_wig`)

WIG files with per-window GC content and mappability scores. **The bin size must match the readcounter window** (pipeline default: 1 Mb). Pre-built hg38 files at 1 Mb resolution (`gc_hg38_1000kb.wig`, `map_hg38_1000kb.wig`) are available from:
https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata

### dbSNP (`--dbsnp`, `--dbsnp_tbi`)

Required by MuSE sump. Must be bgzipped and tabix-indexed:
https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/

### TNScope reference VCFs

| File | Source |
|------|--------|
| `af-only-gnomad.hg38.vcf.gz` | GATK bundle: `gs://gatk-best-practices/somatic-hg38/` |
| `small_exac_common_3.hg38.vcf.gz` | GATK bundle: `gs://gatk-best-practices/somatic-hg38/` |

Both files must be bgzipped and tabix-indexed.

### VEP cache (`--vep_cache`)

```bash
apptainer exec docker://ensemblorg/ensembl-vep:release_114.0 \
    vep_install -a cf -s homo_sapiens -y GRCh38 --CACHE_VERSION 114 -c ~/.vep
```

### Sentieon license and installation

Sentieon reads its license and binary location from environment variables. Set these before running the pipeline (e.g. in your `~/.bashrc` or Slurm job preamble):

```bash
export SENTIEON_LICENSE=/path/to/your.lic
export SENTIEON_INSTALL_DIR=/path/to/sentieon-genomics-202503.xx/bin
export SENTIEON_TMPDIR=/path/to/scratch/sentieon_tmp
```

The pipeline whitelists all three variables so Apptainer passes them into the container (`envWhitelist` in `nextflow.config`). No pipeline parameter is needed.

---

## Output Structure

```
results/
├── bam/
│   ├── tumor/              # Converted tumor BAMs (.bam + .bai)
│   └── normal/             # Converted normal BAMs (.bam + .bai)
├── muse/
│   ├── call/               # MuSE intermediate files (*.MuSE.txt)
│   └── sump/               # MuSE somatic VCFs (*.vcf.gz + .tbi)
├── tnscope/
│   ├── call/               # TNhaplotyper2 raw VCF + auxiliary files
│   │   ├── *_tnhap2-tmp.vcf.gz
│   │   ├── *_orientation
│   │   ├── *_contamination
│   │   └── *_contamination-segments
│   └── filter/             # TNfilter final VCFs (*_tnhap2.vcf.gz + .tbi)
├── vep/                    # Annotated VCFs for both callers (*.vep.vcf.gz + .tbi + .summary.html)
├── hmmcopy/
│   ├── tumor/              # Tumor read count WIG files
│   └── normal/             # Normal read count WIG files
└── ichorcna/
    ├── *.seg               # Copy number segments
    ├── *.cna.seg
    ├── *.seg.txt
    ├── *.params.txt        # Estimated tumor fraction and ploidy
    ├── *.correctedDepth.txt
    ├── *.RData             # R objects for reanalysis
    └── plots/              # Genome-wide CNA plots (PDF)
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

Global limits: 72 CPUs, 512 GB RAM, 168 h per job. Override with a custom config:

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
- **Sentieon** — license and install path set via `SENTIEON_LICENSE` / `SENTIEON_INSTALL_DIR` env vars (see Reference Files)

---

## License

MIT
