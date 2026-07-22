#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Somatic Variant Calling + Copy Number Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CRAM → BAM conversion, MuSE somatic variant calling, and ichorCNA copy number analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT_TUMOR  } from './modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT_NORMAL } from './modules/nf-core/samtools/convert/main'
include { MUSE_CALL                                   } from './modules/nf-core/muse/call/main'
include { MUSE_SUMP                                   } from './modules/nf-core/muse/sump/main'
include { HMMCOPY_READCOUNTER as READCOUNTER_TUMOR    } from './modules/nf-core/hmmcopy/readcounter/main'
include { HMMCOPY_READCOUNTER as READCOUNTER_NORMAL   } from './modules/nf-core/hmmcopy/readcounter/main'
include { ICHORCNA_RUN                                } from './modules/nf-core/ichorcna/run/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_MUSE       } from './modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_VEP_TNSCOPE    } from './modules/nf-core/ensemblvep/vep/main'
include { SENTIEON_TNHAPLOTYPER2                      } from './modules/nf-core/sentieon/tnhaplotyper2/main'
include { SENTIEON_TNFILTER                           } from './modules/nf-core/sentieon/tnfilter/main'
include { GATK4_ASEREADCOUNTER as GATK4_ASEREADCOUNTER_TUMOR  } from './modules/nf-core/gatk4/asereadcounter/main'
include { GATK4_ASEREADCOUNTER as GATK4_ASEREADCOUNTER_NORMAL } from './modules/nf-core/gatk4/asereadcounter/main'
include { HMFTOOLS_AMBER                              } from './modules/nf-core/hmftools/amber/main'
include { HMFTOOLS_COBALT                             } from './modules/nf-core/hmftools/cobalt/main'
include { HMFTOOLS_PURPLE                             } from './modules/nf-core/hmftools/purple/main'

workflow {

    // -------------------------------------------------------------------------
    // Parameter validation
    // -------------------------------------------------------------------------
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // -------------------------------------------------------------------------
    // Reference channels
    // -------------------------------------------------------------------------
    ch_fasta     = channel.fromPath(params.fasta, checkIfExists: true).first()
    ch_fai       = channel.fromPath(params.fai, checkIfExists: true).first()
    ch_dict      = params.dict ? channel.fromPath(params.dict, checkIfExists: true).first() : channel.value([])
    ch_dbsnp     = channel.fromPath(params.dbsnp, checkIfExists: true).first()
    ch_dbsnp_tbi = channel.fromPath(params.dbsnp_tbi, checkIfExists: true).first()
    ch_vep_cache          = params.vep_cache          ? channel.fromPath(params.vep_cache, checkIfExists: true).first()          : channel.value([])
    ch_het_pon            = params.het_pon            ? channel.fromPath(params.het_pon, checkIfExists: true).first()            : channel.value([])
    ch_het_pon_tbi        = params.het_pon_tbi        ? channel.fromPath(params.het_pon_tbi, checkIfExists: true).first()        : channel.value([])
    ch_gc_profile         = params.gc_profile         ? channel.fromPath(params.gc_profile, checkIfExists: true).first()         : channel.value([])
    ch_ensembl_data_dir   = params.ensembl_data_dir   ? channel.fromPath(params.ensembl_data_dir, checkIfExists: true).first()   : channel.value([])
    ch_germline_vcf       = params.germline_vcf       ? channel.fromPath(params.germline_vcf, checkIfExists: true).first()       : channel.value([])
    ch_germline_vcf_tbi   = params.germline_vcf_tbi   ? channel.fromPath(params.germline_vcf_tbi, checkIfExists: true).first()   : channel.value([])
    ch_contamination_vcf  = params.contamination_vcf  ? channel.fromPath(params.contamination_vcf, checkIfExists: true).first()  : channel.value([])
    ch_contamination_vcf_tbi = params.contamination_vcf_tbi ? channel.fromPath(params.contamination_vcf_tbi, checkIfExists: true).first() : channel.value([])
    ch_gc_wig    = params.gc_wig         ? channel.fromPath(params.gc_wig, checkIfExists: true).first()         : channel.value([])
    ch_map_wig   = params.map_wig        ? channel.fromPath(params.map_wig, checkIfExists: true).first()        : channel.value([])
    ch_centromere = params.centromere     ? channel.fromPath(params.centromere, checkIfExists: true).first()     : channel.value([])
    ch_pon       = params.panel_of_normals ? channel.fromPath(params.panel_of_normals, checkIfExists: true).first() : channel.value([])
    ch_rep_time  = params.rep_time_wig   ? channel.fromPath(params.rep_time_wig, checkIfExists: true).first()   : channel.value([])
    ch_exons     = params.exons          ? channel.fromPath(params.exons, checkIfExists: true).first()          : channel.value([])

    // Reference tuple for samtools/convert: [meta2, fasta, fai]
    ch_ref = ch_fasta.combine(ch_fai).map { fasta, fai ->
        tuple([id: 'reference'], fasta, fai)
    }

    // Reference tuple for hmmcopy/readcounter: [meta2, fasta]
    ch_ref_fasta = ch_fasta.map { fasta -> tuple([id: 'reference'], fasta) }

    // -------------------------------------------------------------------------
    // Parse samplesheet
    // Required: patient, tumor_cram, tumor_crai, normal_cram, normal_crai
    // Optional: tumor_sample, normal_sample (SM tags in CRAM read groups;
    //           required for TNScope — must match @RG SM: in the CRAM header)
    // -------------------------------------------------------------------------
    ch_input = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id:            row.patient,
                // Strip leading "SM:" if the user copied it from samtools @RG output
                tumor_sample:  (row.tumor_sample  ?: row.patient).replaceAll(/^SM:/, ''),
                normal_sample: (row.normal_sample ?: (row.patient + '_NORMAL')).replaceAll(/^SM:/, ''),
            ]
            tuple(
                meta,
                file(row.tumor_cram, checkIfExists: true),
                file(row.tumor_crai, checkIfExists: true),
                file(row.normal_cram, checkIfExists: true),
                file(row.normal_crai, checkIfExists: true)
            )
        }

    // Split into tumor and normal channels.
    // Normal gets a distinct meta.id (_N suffix) so samtools/convert and
    // hmmcopy/readcounter output different filenames from the tumor. The
    // original patient id is preserved in meta.patient for re-joining later.
    ch_tumor_cram = ch_input.map { meta, tumor_cram, tumor_crai, _normal_cram, _normal_crai ->
        tuple(meta, tumor_cram, tumor_crai)
    }

    ch_normal_cram = ch_input.map { meta, _tumor_cram, _tumor_crai, normal_cram, normal_crai ->
        tuple([id: "${meta.id}_N", patient: meta.id], normal_cram, normal_crai)
    }

    // -------------------------------------------------------------------------
    // Step 1: CRAM → BAM conversion
    // -------------------------------------------------------------------------

    /*
    Convert tumor CRAMs to BAM format for downstream tools that require BAM input.
    samtools/convert input: tuple(meta, input, index), tuple(meta2, fasta, fai)
    */
    SAMTOOLS_CONVERT_TUMOR(ch_tumor_cram, ch_ref)

    /*
    Convert matched normal CRAMs to BAM format.
    */
    SAMTOOLS_CONVERT_NORMAL(ch_normal_cram, ch_ref)

    // Collect tumor BAM + BAI (full meta — tumor_sample/normal_sample needed by AMBER/COBALT)
    ch_tumor_bam = SAMTOOLS_CONVERT_TUMOR.out.bam
        .join(SAMTOOLS_CONVERT_TUMOR.out.bai)

    // Collect normal BAM + BAI; remap meta back to patient id for downstream joins
    ch_normal_bam = SAMTOOLS_CONVERT_NORMAL.out.bam
        .join(SAMTOOLS_CONVERT_NORMAL.out.bai)
        .map { meta, bam, bai -> tuple([id: meta.patient], bam, bai) }

    // Joined tumor+normal BAM channel with full tumor meta preserved.
    // Join on bare string id so the different-shaped meta maps don't block the match.
    ch_tumor_normal_bam = ch_tumor_bam
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(
            SAMTOOLS_CONVERT_NORMAL.out.bam
                .join(SAMTOOLS_CONVERT_NORMAL.out.bai)
                .map { meta, bam, bai -> tuple(meta.patient, bam, bai) }
        )
        .map { _id, full_meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(full_meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }

    // id-only meta version for ASE tumor (joins with id-only VCF channel)
    ch_tumor_bam_id = ch_tumor_bam
        .map { meta, bam, bai -> tuple([id: meta.id], bam, bai) }

    // -------------------------------------------------------------------------
    // Step 2: MuSE somatic variant calling
    // -------------------------------------------------------------------------

    /*
    Combine tumor and normal BAMs with the reference FASTA for MuSE call.
    muse/call input: tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, reference)
    */
    ch_muse_input = ch_tumor_normal_bam
        .combine(ch_fasta)
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, fasta ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, fasta)
        }

    MUSE_CALL(ch_muse_input)

    /*
    Run MuSE sump to compute tier-based cutoffs and finalize somatic variants.
    Mode is controlled by params.wgs via ext.args: -G (WGS) or -E (WXS).
    muse/sump input: tuple(meta, muse_call_txt, ref_vcf, ref_vcf_tbi)
    */
    ch_sump_input = MUSE_CALL.out.txt
        .combine(ch_dbsnp)
        .combine(ch_dbsnp_tbi)
        .map { meta, txt, dbsnp, dbsnp_tbi ->
            tuple(meta, txt, dbsnp, dbsnp_tbi)
        }

    MUSE_SUMP(ch_sump_input)

    // -------------------------------------------------------------------------
    // Step 3: ichorCNA copy number analysis
    // -------------------------------------------------------------------------

    /*
    Count reads in genomic windows for the tumor BAM.
    hmmcopy/readcounter input: tuple(meta, bam, bai), tuple(meta2, fasta)
    */
    READCOUNTER_TUMOR(ch_tumor_bam, ch_ref_fasta)

    /*
    Count reads in genomic windows for the matched normal BAM.
    */
    READCOUNTER_NORMAL(ch_normal_bam, ch_ref_fasta)

    /*
    Run ichorCNA using the tumor wig, normal wig, GC content, and mappability data
    to estimate copy number alterations and tumor fraction.
    ichorcna/run input: tuple(meta, wig), gc_wig, map_wig, normal_wig, normal_background, centromere, rep_time_wig, exons
    */
    // Join normal wig to tumor wig by patient id so each tumor sample gets its
    // paired normal wig rather than relying on channel ordering.
    ch_ichorcna_input = READCOUNTER_TUMOR.out.wig
        .map { meta, wig -> tuple(meta.id, meta, wig) }
        .join(
            READCOUNTER_NORMAL.out.wig
                .map { meta, wig -> tuple(meta.patient, wig) }
        )
        .map { _id, meta, tumor_wig, normal_wig -> tuple(meta, tumor_wig, normal_wig) }

    ch_ichorcna_input
        .multiMap { meta, tumor_wig, normal_wig ->
            tumor:  tuple(meta, tumor_wig)
            normal: normal_wig
        }
        .set { ch_ichorcna_split }

    ICHORCNA_RUN(
        ch_ichorcna_split.tumor,
        ch_gc_wig,
        ch_map_wig,
        ch_ichorcna_split.normal,
        ch_pon,
        ch_centromere,
        ch_rep_time,
        ch_exons
    )

    // -------------------------------------------------------------------------
    // Step 4: TNScope somatic variant calling (Sentieon)
    // -------------------------------------------------------------------------

    /*
    Run TNhaplotyper2 + OrientationBias + ContaminationModel in a single
    sentieon driver call. Operates on CRAMs directly (no BAM conversion needed).
    meta.tumor_sample and meta.normal_sample must match the SM tags in the
    CRAM read groups (add tumor_sample/normal_sample columns to the samplesheet).
    sentieon/tnhaplotyper2 input: tuple(meta, tumor_cram, tumor_crai, normal_cram, normal_crai),
        fasta, fai, germline_vcf, germline_vcf_tbi, contamination_vcf, contamination_vcf_tbi, license
    */
    ch_tnscope_input = ch_input

    SENTIEON_TNHAPLOTYPER2(
        ch_tnscope_input,
        ch_fasta,
        ch_fai,
        ch_germline_vcf,
        ch_germline_vcf_tbi,
        ch_contamination_vcf,
        ch_contamination_vcf_tbi
    )

    /*
    Filter the raw TNhaplotyper2 VCF using orientation bias and contamination
    estimates to produce the final somatic VCF.
    */
    ch_tnfilter_input = SENTIEON_TNHAPLOTYPER2.out.vcf_tmp
        .join(SENTIEON_TNHAPLOTYPER2.out.vcf_tmp_tbi)
        .join(SENTIEON_TNHAPLOTYPER2.out.vcf_tmp_stats)
        .join(SENTIEON_TNHAPLOTYPER2.out.orientation)
        .join(SENTIEON_TNHAPLOTYPER2.out.contamination)
        .join(SENTIEON_TNHAPLOTYPER2.out.contamination_segments)

    SENTIEON_TNFILTER(
        ch_tnfilter_input,
        ch_fasta,
        ch_fai
    )

    // -------------------------------------------------------------------------
    // Step 5: VEP annotation — MuSE and TNScope VCFs
    // -------------------------------------------------------------------------

    /*
    Annotate somatic VCFs from both callers with Ensembl VEP using a local
    offline cache. The two VCF channels are mixed so VEP runs once per sample
    per caller. Skipped when --vep_cache is not supplied.
    ensemblvep/vep input: tuple(meta, vcf, tbi), cache_dir, genome, cache_version
    */
    ch_muse_for_vep    = MUSE_SUMP.out.vcf.join(MUSE_SUMP.out.tbi)
    ch_tnscope_for_vep = SENTIEON_TNFILTER.out.vcf.join(SENTIEON_TNFILTER.out.tbi)

    ENSEMBLVEP_VEP_MUSE(
        ch_muse_for_vep,
        ch_vep_cache,
        params.vep_genome,
        params.vep_cache_version
    )

    ENSEMBLVEP_VEP_TNSCOPE(
        ch_tnscope_for_vep,
        ch_vep_cache,
        params.vep_genome,
        params.vep_cache_version
    )

    // -------------------------------------------------------------------------
    // Step 6: Purple — tumor purity and ploidy estimation
    // -------------------------------------------------------------------------

    /*
    AMBER computes B-allele frequencies from the tumor and matched normal BAMs.
    hmftools/amber input: tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai),
        fasta, fai, het_pon, het_pon_tbi, ref_genome_version
    */
    ch_purple_bam_input = ch_tumor_normal_bam

    HMFTOOLS_AMBER(
        ch_purple_bam_input,
        ch_fasta,
        ch_fai,
        ch_het_pon,
        ch_het_pon_tbi,
        params.ref_genome_version
    )

    /*
    COBALT computes read-depth ratios from the tumor and matched normal BAMs.
    hmftools/cobalt input: tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai),
        fasta, fai, gc_profile
    */
    HMFTOOLS_COBALT(
        ch_purple_bam_input,
        ch_fasta,
        ch_fai,
        ch_gc_profile,
        params.ref_genome_version
    )

    /*
    PURPLE estimates tumor purity and ploidy from AMBER BAFs, COBALT read
    ratios, and the TNScope somatic VCF. The somatic VCF improves SNV
    copy-number context annotation but is optional.
    hmftools/purple input: tuple(meta, amber_dir, cobalt_dir, somatic_vcf, somatic_vcf_tbi),
        fasta, fai, gc_profile, ref_genome_version
    */
    ch_purple_input = HMFTOOLS_AMBER.out.amber_dir
        .map { meta, dir -> tuple(meta.id, meta, dir) }
        .join(HMFTOOLS_COBALT.out.cobalt_dir.map { meta, dir -> tuple(meta.id, dir) })
        .join(SENTIEON_TNFILTER.out.vcf.map      { meta, vcf -> tuple(meta.id, vcf) })
        .join(SENTIEON_TNFILTER.out.tbi.map      { meta, tbi -> tuple(meta.id, tbi) })
        .map { _id, meta, amber_dir, cobalt_dir, vcf, tbi -> tuple(meta, amber_dir, cobalt_dir, vcf, tbi) }

    HMFTOOLS_PURPLE(
        ch_purple_input,
        ch_fasta,
        ch_fai,
        ch_gc_profile,
        params.ref_genome_version,
        ch_ensembl_data_dir
    )

    // -------------------------------------------------------------------------
    // Step 7: ASEReadCounter — allele-specific read counts at somatic SNV sites
    // -------------------------------------------------------------------------

    /*
    SelectVariants filters the TNScope VCF to passing biallelic SNPs, then
    ASEReadCounter measures the allelic read depth at each site in the tumor
    and matched normal BAMs. The meta2 tuple carries the VCF so it is
    broadcast to every BAM sample by key.
    gatk4/asereadcounter input: tuple(meta, bam, bai), tuple(meta2, vcf, vcf_tbi), fasta, fai
    */
    // Normalize TNScope VCF meta to id-only so it joins with BAM channels
    ch_ase_vcf = SENTIEON_TNFILTER.out.vcf
        .join(SENTIEON_TNFILTER.out.tbi)
        .map { meta, vcf, tbi -> tuple([id: meta.id], vcf, tbi) }

    ch_ase_tumor_input = ch_tumor_bam_id
        .join(ch_ase_vcf)
        .map { meta, bam, bai, vcf, tbi -> tuple(meta, bam, bai, tuple(meta, vcf, tbi)) }

    GATK4_ASEREADCOUNTER_TUMOR(
        ch_ase_tumor_input.map { meta, bam, bai, _vcf_tuple -> tuple(meta, bam, bai) },
        ch_ase_tumor_input.map { _meta, _bam, _bai, vcf_tuple -> vcf_tuple },
        ch_fasta,
        ch_fai,
        ch_dict
    )

    ch_ase_normal_input = ch_normal_bam
        .map { meta, bam, bai -> tuple(meta.id, bam, bai) }
        .join(ch_ase_vcf.map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi) })
        .map { id, bam, bai, vcf, tbi -> tuple([id: id], bam, bai, tuple([id: id], vcf, tbi)) }

    GATK4_ASEREADCOUNTER_NORMAL(
        ch_ase_normal_input.map { meta, bam, bai, _vcf_tuple -> tuple(meta, bam, bai) },
        ch_ase_normal_input.map { _meta, _bam, _bai, vcf_tuple -> vcf_tuple },
        ch_fasta,
        ch_fai,
        ch_dict
    )
}
