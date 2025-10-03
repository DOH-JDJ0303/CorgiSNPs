//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTERQDUMP  } from '../../../modules/local/fasterq-dump/main'
include { SEQTK_SAMPLE } from '../../../modules/local/seqtk/sample/main'
include { FASTQC       } from '../../../modules/nf-core/fastqc/main'
include { FASTP        } from '../../../modules/nf-core/fastp/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE {

    take:
    // Channel: [ val(meta), path(reads) ] plus optional fields (sp, sb, ref, sra)
    ch_samplesheet

    main:

    // Collectors
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // -------------------------------------------------------------------------
    // Branch: separate SRA vs non-SRA rows from samplesheet
    // Produces ch_sra_branch.sra_true and ch_sra_branch.sra_false
    // -------------------------------------------------------------------------
    ch_samplesheet
        .branch { meta, reads, sp, sb, ref, sra ->
            sra_true :  sra
            sra_false: !sra
        }
        .set { ch_sra_branch }

    // -------------------------------------------------------------------------
    // MODULE: Download reads from SRA for sra_true rows
    // Input reshaped to [ meta, sra ]
    // -------------------------------------------------------------------------
    FASTERQDUMP(
        ch_sra_branch.sra_true.map { meta, reads, sp, sb, ref, sra -> [ meta, sra ] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)

    // -------------------------------------------------------------------------
    // Merge SRA-derived reads back with pass-through non-SRA reads
    // Ensure meta.single_end is set based on number of read files
    // Output ch_samplesheet retains shape: [ meta, reads, sp, sb, ref ]
    // -------------------------------------------------------------------------
    ch_sra_branch
        .sra_true
        .join(FASTERQDUMP.out.reads, by: 0)
        .map   { meta, empty_reads, sp, sb, ref, sra, reads -> [ meta, reads, sp, sb, ref ] }
        .concat( ch_sra_branch.sra_false.map { meta, reads, sp, sb, ref, sra -> [ meta, reads, sp, sb, ref ] } )
        .map { meta, reads, sp, sb, ref ->
            meta.single_end = (reads.size() == 1)
            [ meta, reads, sp, sb, ref ]
        }
        .set { ch_samplesheet }

    // Split into read and meta streams for downstream modules
    ch_samplesheet.map { meta, reads, sp, sb, ref -> [ meta, reads ] }.set { ch_reads }
    ch_samplesheet.map { meta, reads, sp, sb, ref -> [ meta, sp, sb, ref ] }.set { ch_meta }

    // -------------------------------------------------------------------------
    // MODULE: Downsample reads with seqtk sample (if --max_reads provided)
    // Keeps tuple shape [ meta, reads ] throughout.
    // -------------------------------------------------------------------------
    if (params.max_reads) {

        // Compute total read count (approx) by counting first R1 and doubling (paired)
        ch_reads
            .map { meta, reads -> [ meta, reads, reads[0].countFastq() * (meta.single_end ? 1 : 2) ] }
            .branch { meta, reads, n ->
                ok  : n <= params.max_reads
                high: n >  params.max_reads
            }
            .set { ch_reads }

        // For high-coverage samples, sample each mate independently
        SEQTK_SAMPLE(
            ch_reads
                .high
                .transpose()                                  // [ meta, read, n ]
                .map { meta, read, n -> [ meta, read, params.max_reads ] }
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

        // Re-assemble paired reads and merge with ok set
        SEQTK_SAMPLE.out.read
            .groupTuple(by: 0)                                 // -> [ meta, [paths...] ]
            .concat( ch_reads.ok.map { meta, reads, n -> [ meta, reads ] } )
            .set { ch_reads }
    }

    // -------------------------------------------------------------------------
    // MODULE: FastQC (adds zips to MultiQC input and versions to collector)
    // -------------------------------------------------------------------------
    FASTQC(
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect { it[1] } )
    ch_versions      = ch_versions     .mix( FASTQC.out.versions.first() )

    // -------------------------------------------------------------------------
    // MODULE: Fastp (trimming/filters + JSON stats). Replace reads with trimmed.
    // -------------------------------------------------------------------------
    FASTP(
        ch_reads,
        [],
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    FASTP.out.reads.set { ch_reads }

    emit:
    meta          = ch_meta
    reads         = ch_reads
    read_stats    = FASTP.out.json
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
