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
    // MODULE: Download reads from SRA for sra_true rows
    // Input reshaped to [ meta, sra ]
    // -------------------------------------------------------------------------
    FASTERQDUMP(
        ch_samplesheet
            .filter{ it.meta.sra }
            .map { it -> [ it.meta, it.meta.sra ] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)

    // -------------------------------------------------------------------------
    // Merge SRA-derived reads back with pass-through non-SRA reads
    // Ensure meta.single_end is set based on number of read files
    // Output ch_samplesheet retains shape: [ meta, reads, sp, sb, ref ]
    // -------------------------------------------------------------------------
    FASTERQDUMP
        .out
        .reads
        .map{ meta, reads -> 
            def new_meta = meta + [single_end: reads.size() == 1]
            [meta: new_meta, reads: reads] 
        }
        .concat( 
            ch_samplesheet.filter{ ! it.meta.sra }
        )
        .map{ it ->
            def new_meta = it.meta.findAll{ k, v -> k != 'sra' }
            return it + [meta: new_meta]
         }
        .set { ch_samplesheet }

    // -------------------------------------------------------------------------
    // MODULE: Downsample reads with seqtk sample (if --max_reads provided)
    // Keeps tuple shape [ meta, reads ] throughout.
    // -------------------------------------------------------------------------
    if (params.max_reads) {

        // Compute total read count (approx) by counting first R1 and doubling (paired)
        ch_samplesheet
            .map { it -> [it, it.reads[0].countFastq() * (it.meta.single_end ? 1 : 2) ] }
            .branch { it, n ->
                ok  : n <= params.max_reads
                high: n >  params.max_reads
            }
            .set { ch_samplesheet }

        // For high-coverage samples, sample each mate independently
        SEQTK_SAMPLE(
            ch_samplesheet
                .high
                .map{ it, n -> [it.meta, it.reads, params.max_reads] }
                .transpose()
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

        // Re-assemble paired reads and merge with ok set
        SEQTK_SAMPLE
            .out
            .read
            .groupTuple(by: 0)
            .map{ meta, reads -> [meta: meta, reads: reads] }                            // -> [ meta, [paths...] ]
            .concat( ch_samplesheet.ok.map { it, n -> it } )
            .set { ch_samplesheet }
    }

    ch_samplesheet = ch_samplesheet.map{ [it.meta, it.reads] }

    // -------------------------------------------------------------------------
    // MODULE: FastQC (adds zips to MultiQC input and versions to collector)
    // -------------------------------------------------------------------------
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect { it[1] } )
    ch_versions      = ch_versions     .mix( FASTQC.out.versions.first() )

    // -------------------------------------------------------------------------
    // MODULE: Fastp (trimming/filters + JSON stats). Replace reads with trimmed.
    // -------------------------------------------------------------------------
    FASTP(
        ch_samplesheet,
        [],
        false,
        false,
        false
    )
    ch_versions    = ch_versions.mix(FASTP.out.versions.first())
    ch_samplesheet = FASTP.out.reads

    emit:
    samplesheet   = ch_samplesheet
    read_stats    = FASTP.out.json
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
