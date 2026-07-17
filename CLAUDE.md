# CLAUDE.md — moonshot_pipeline

## What this pipeline does

Somatic variant calling and copy number analysis from matched tumor/normal CRAM pairs.

**Workflow (`main.nf`):**
1. **CRAM → BAM** — `samtools/convert` for tumor and normal independently
2. **Somatic SNVs (MuSE)** — `muse/call` then `muse/sump` (WGS `-G` or WXS `-E` mode)
3. **Copy number (ichorCNA)** — `hmmcopy/readcounter` on both BAMs → `ichorcna/run`
4. **Somatic SNVs (TNScope)** — `sentieon/tnhaplotyper2` then `sentieon/tnfilter` (requires Sentieon license)
5. **VEP annotation** — `ensemblvep/vep` annotates both MuSE and TNScope VCFs (requires `--vep_cache`)
6. **Tumor purity (Purple)** — `hmftools/amber` → `hmftools/cobalt` → `hmftools/purple`
7. **Allele-specific expression** — `gatk4/asereadcounter` counts reads at biallelic somatic SNP sites in both tumor and normal BAMs

## Running the pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --fasta /path/to/ref.fa \
  --fai   /path/to/ref.fa.fai \
  --dbsnp /path/to/dbsnp.vcf.gz \
  --dbsnp_tbi /path/to/dbsnp.vcf.gz.tbi \
  --gc_wig  /path/to/gc.wig \
  --map_wig /path/to/map.wig \
  --outdir results
```

Default executor is **Slurm + Apptainer**. Use `-profile local` to run locally.
Use `-profile wxs` for whole-exome mode (switches MuSE sump to `-E`).

## Samplesheet format (`assets/samplesheet.csv`)

```
patient,tumor_cram,tumor_crai,normal_cram,normal_crai
PATIENT_01,/path/tumor.cram,/path/tumor.cram.crai,/path/normal.cram,/path/normal.cram.crai
```

One row per patient. All paths must be absolute or resolvable from the run directory.

## Key parameters (`nextflow.config`)

| Parameter | Required | Description |
|---|---|---|
| `--input` | yes | Samplesheet CSV |
| `--fasta` / `--fai` | yes | Reference genome + index |
| `--dict` | no | Reference sequence dictionary (.dict) — required by ASEReadCounter |
| `--dbsnp` / `--dbsnp_tbi` | yes | dbSNP VCF (bgzipped + tabix) |
| `--gc_wig` | yes* | GC content wig (hmmcopy gcCounter) |
| `--map_wig` | yes* | Mappability wig (hmmcopy mapCounter) |
| `--centromere` | no | Centromere locations for ichorCNA |
| `--panel_of_normals` | no | ichorCNA panel of normals RDS |
| `--rep_time_wig` | no | Replication timing wig |
| `--exons` | no | Exon BED for annotation |
| `--wgs` | no | `true` (default) = WGS, `false` = WXS |
| `--het_pon` | no | HMF GermlineHetPon VCF for AMBER |
| `--het_pon_tbi` | no | Tabix index for het PON |
| `--gc_profile` | no | HMF GC profile CNP for COBALT + PURPLE |
| `--ref_genome_version` | no | HMF ref genome version (default `V38`) |
| `--vep_cache` | no | Path to VEP cache directory (`~/.vep`) |
| `--vep_genome` | no | Assembly name for VEP (default `GRCh38`) |
| `--vep_cache_version` | no | VEP cache version (default `114`) |
| `--outdir` | no | Output directory (default `./results`) |

*Required for ichorCNA to produce meaningful output.

## Module structure

All modules are nf-core style under `modules/nf-core/`:

```
modules/nf-core/
  samtools/convert/      # CRAM → BAM
  muse/call/             # MuSE somatic calling
  muse/sump/             # MuSE tier cutoffs → VCF
  ensemblvep/vep/        # VEP annotation of somatic VCFs
  hmmcopy/readcounter/   # Read counting → WIG
  ichorcna/run/          # CNA + tumor fraction estimation
  sentieon/tnhaplotyper2/  # TNhaplotyper2 + OrientationBias + ContaminationModel
  sentieon/tnfilter/       # TNfilter → final somatic VCF
  hmftools/amber/          # B-allele frequencies (Purple prereq)
  hmftools/cobalt/         # Read-depth ratios (Purple prereq)
  hmftools/purple/         # Tumor purity, ploidy, copy number
  gatk4/asereadcounter/    # Allele-specific read counts at somatic SNP sites
```

Each module has `main.nf`, `meta.yml`, `environment.yml`, and `tests/`.

## Output structure

```
results/
  bam/tumor/             # Tumor BAMs
  bam/normal/            # Normal BAMs
  muse/call/             # *.MuSE.txt intermediate files
  muse/sump/             # *.vcf.gz somatic SNV calls
  muse/annotation/       # VEP-annotated MuSE VCFs
  tnscope/annotation/    # VEP-annotated TNScope VCFs
  hmmcopy/tumor/         # Tumor read count WIG
  hmmcopy/normal/        # Normal read count WIG
  ase/tumor/             # *.ase.csv allele counts (tumor BAM)
  ase/normal/            # *.ase.csv allele counts (normal BAM)
  ichorcna/              # Segments, params, corrected depth, RData
  ichorcna/plots/        # PDF plots
```

## HPC / execution

- **Executor**: Slurm (default), overridable with `-profile local`
- **Containers**: Apptainer (Singularity), cached at `~/.apptainer/cache`; pulled via Wave from conda environments
- **Resource labels**: `process_single` → `process_high` (1–12 CPUs, 1–64 GB RAM)
- **Limits**: 72 CPUs, 512 GB RAM, 168 h per job
- **Read counter window**: 1 Mb (`-w 1000000`), MAPQ ≥ 20 (`-q 20`)

## Plugin

Uses `nf-schema@2.7.2` for parameter validation (`validateParameters()`) and summary logging.
Requires Nextflow ≥ 24.04.0.

## Making changes

- **Add a new tool**: drop an nf-core-style module under `modules/nf-core/<tool>/`, add `include` + call in `main.nf`, add `withName` block in `conf/modules.config`.
- **Change MuSE mode**: set `--wgs false` (WXS) or `--wgs true` (WGS); the `ext.args` in `modules.config` handles the flag automatically.
- **Change read counter resolution**: edit `ext.args` for `READCOUNTER_TUMOR`/`READCOUNTER_NORMAL` in `conf/modules.config` (`-w` flag).
- **Test locally**: `-profile local` disables Slurm; containers still pulled via Wave/Apptainer unless you also override the container engine.
