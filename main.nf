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
include { ENSEMBLVEP_VEP                              } from './modules/nf-core/ensemblvep/vep/main'
include { SENTIEON_TNHAPLOTYPER2                      } from './modules/nf-core/sentieon/tnhaplotyper2/main'
include { SENTIEON_TNFILTER                           } from './modules/nf-core/sentieon/tnfilter/main'

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
    ch_dbsnp     = channel.fromPath(params.dbsnp, checkIfExists: true).first()
    ch_dbsnp_tbi = channel.fromPath(params.dbsnp_tbi, checkIfExists: true).first()
    ch_vep_cache          = params.vep_cache          ? channel.fromPath(params.vep_cache, checkIfExists: true).first()          : channel.value([])
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
                tumor_sample:  row.tumor_sample  ?: row.patient,
                normal_sample: row.normal_sample ?: (row.patient + '_NORMAL'),
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

    // Collect tumor BAM + BAI
    ch_tumor_bam = SAMTOOLS_CONVERT_TUMOR.out.bam
        .join(SAMTOOLS_CONVERT_TUMOR.out.bai)

    // Collect normal BAM + BAI; remap meta back to patient id for downstream joins
    ch_normal_bam = SAMTOOLS_CONVERT_NORMAL.out.bam
        .join(SAMTOOLS_CONVERT_NORMAL.out.bai)
        .map { meta, bam, bai -> tuple([id: meta.patient], bam, bai) }

    // -------------------------------------------------------------------------
    // Step 2: MuSE somatic variant calling
    // -------------------------------------------------------------------------

    /*
    Combine tumor and normal BAMs with the reference FASTA for MuSE call.
    muse/call input: tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, reference)
    */
    ch_muse_input = ch_tumor_bam
        .join(ch_normal_bam)
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
    ch_normal_wig_keyed = READCOUNTER_NORMAL.out.wig
        .map { meta, wig -> tuple([id: meta.patient], wig) }

    ch_ichorcna_input = READCOUNTER_TUMOR.out.wig
        .join(ch_normal_wig_keyed)
        .map { meta, tumor_wig, normal_wig -> tuple(meta, tumor_wig, normal_wig) }

    ICHORCNA_RUN(
        ch_ichorcna_input.map { meta, tumor_wig, _normal_wig -> tuple(meta, tumor_wig) },
        ch_gc_wig,
        ch_map_wig,
        ch_ichorcna_input.map { _meta, _tumor_wig, normal_wig -> normal_wig },
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
        .map { meta, vcf, tbi -> tuple(meta + [caller: 'muse'], vcf, tbi) }

    ch_tnscope_for_vep = SENTIEON_TNFILTER.out.vcf.join(SENTIEON_TNFILTER.out.tbi)
        .map { meta, vcf, tbi -> tuple(meta + [caller: 'tnscope'], vcf, tbi) }

    ch_vep_input = ch_muse_for_vep.mix(ch_tnscope_for_vep)

    ENSEMBLVEP_VEP(
        ch_vep_input,
        ch_vep_cache,
        params.vep_genome,
        params.vep_cache_version
    )
}
