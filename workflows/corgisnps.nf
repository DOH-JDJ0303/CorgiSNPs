import Utils

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_corgisnps_pipeline'

// Subworkflows
include { PREPARE  } from '../subworkflows/local/prepare'
include { CLASSIFY } from '../subworkflows/local/classify'
include { VARIANTS } from '../subworkflows/local/variants'
include { AMR      } from '../subworkflows/local/amr'
include { PHYLO    } from '../subworkflows/local/phylo'

// Report / Summary Modules
include { SUMMARYLINE    } from '../modules/local/report/main'
include { ADD_AMR        } from '../modules/local/report/main'
include { REPORT_SPECIES } from '../modules/local/report/main'
include { REPORT_ALL     } from '../modules/local/report/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CORGISNPS {

    take:
    // Channel: samplesheet read in from --input
    ch_samplesheet

    main:

    // Collectors for versions and MultiQC input files
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ---------------------------
    // PREPARE
    // ---------------------------
    PREPARE(
        ch_samplesheet
    )
    ch_versions       = ch_versions.mix(PREPARE.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(PREPARE.out.multiqc_files)
    PREPARE.out.reads     .set{ ch_reads }
    PREPARE.out.meta      .set{ ch_meta  }
    PREPARE.out.read_stats.set{ ch_read_stats }

    // Initialize empty downstream channels.
    ch_blank     = ch_reads.map{ [ it[0], [] ] }
    ch_denovo    = ch_blank
    ch_species   = ch_blank
    ch_subtype   = ch_blank
    ch_depth     = ch_blank
    ch_amr       = ch_blank
    ch_aln_stats = ch_blank
    ch_tree      = ch_blank
    ch_dist      = ch_blank

    // ---------------------------
    // CLASSIFY (optional)
    // ---------------------------
    if (params.classify) {

        // Branch meta tuples into those needing classification and those already classified.
        // classify    : items where either species (sp) or subtype (sb) is missing
        // no_classify : items with both sp and sb present
        ch_meta
            .branch { meta, sp, sb, ref ->
                classify   : !(sp && sb)
                no_classify:  (sp && sb)
            }
            .set { ch_meta_branch }

        CLASSIFY(
            // Join reads with only the meta needing classification (by matching meta key)
            ch_reads.join( ch_meta_branch.classify.map{ [ it[0] ] } ),
            ch_meta_branch.classify
        )
        ch_versions = ch_versions.mix(CLASSIFY.out.versions)

        // Merge back classified meta with pass-through meta
        CLASSIFY.out.meta.concat(ch_meta_branch.no_classify).set{ ch_meta }
        CLASSIFY.out.denovo  .set{ ch_denovo }
        CLASSIFY.out.species .set{ ch_species }
        CLASSIFY.out.subtype .set{ ch_subtype }
    }

    // Sanitize species/subtype strings in-flight
    ch_meta.map { meta, sp, sb, ref -> [ meta, Utils.sanitize(sp), Utils.sanitize(sb), ref ] }.set{ch_meta}


    // -------------------------------------------------------------------------
    // Collate per-sample inputs for SUMMARYLINE.
    // Joins preserve samples lacking some inputs via 'remainder: true';
    // missing items are replaced with [] to keep tuple shapes consistent.
    // -------------------------------------------------------------------------
    ch_meta.map{meta, sp, sb, ref -> [meta]}
        .join(ch_read_stats, remainder: true)
        .join(ch_denovo,     remainder: true)
        .join(ch_species,    remainder: true)
        .join(ch_subtype,    remainder: true)
        .map { it.collect { k -> k ? k : [] } }
        .set { ch_samples }

    // -------------------------------------------------------------------------
    // Per-sample summary lines
    // -------------------------------------------------------------------------
    SUMMARYLINE(
        ch_samples,
        file(params.input),
        file(params.ncbi_stats)
    )
    ch_versions = ch_versions.mix(SUMMARYLINE.out.versions.first())

    SUMMARYLINE
        .out
        .summary
        .splitCsv(header: true)
        .map{ meta, data -> [meta, params.auto_qc ? (data.containsKey('qc_status') ? data['qc_status'] == 'PASS' : false) : true ] }
        .branch{ meta, status ->
            pass: status
            not_pass: !status }
        .set{ch_auto_qc}

    ch_reads.join(ch_auto_qc.pass.map{it[0]}).set{ch_reads_pass}
    ch_meta.join(ch_auto_qc.pass.map{it[0]}).set{ch_meta_pass}

    SUMMARYLINE.out.summary.join(ch_auto_qc.pass.map{it[0]}).set{ch_summary_pass}
    SUMMARYLINE.out.summary.join(ch_auto_qc.not_pass.map{it[0]}).set{ch_summary_not_pass}
    // ---------------------------
    // VARIANTS / AMR / PHYLO (optional)
    // ---------------------------
    if (params.variants) {
        VARIANTS(
            ch_reads_pass,
            ch_meta_pass,
            true
        )
        ch_versions = ch_versions.mix(VARIANTS.out.versions)
        VARIANTS.out.meta .set{ ch_meta_pass }
        VARIANTS.out.depth.set{ ch_depth }

        if(params.amr){
            AMR(
                VARIANTS.out.meta,
                VARIANTS.out.bam,
                VARIANTS.out.vcf
            )
            ch_versions = ch_versions.mix(AMR.out.versions)
            
            ch_summary_pass
                .join(AMR.out.summary, remainder: true)
                .branch{ meta, summaryline, amr_summary -> 
                    db_exists: amr_summary
                    db_miss: !amr_summary  }
                .set{ ch_summary_pass_amr }

            ADD_AMR(
                ch_summary_pass_amr.db_exists
            )
            ch_versions = ch_versions.mix(ADD_AMR.out.versions)
            ADD_AMR.out.summary
                .concat(ch_summary_pass_amr.db_miss.map{[it[0],it[1]]})
                .set{ ch_summary_pass }
        }
        if(params.phylo){
            PHYLO(
                VARIANTS.out.aln,
                ch_meta_pass
            )
            ch_versions = ch_versions.mix(PHYLO.out.versions)
            PHYLO.out.aln_stats.set{ ch_aln_stats }
            PHYLO.out.tree     .set{ ch_tree }
            PHYLO.out.dist     .set{ ch_dist }
        }
    }
    
    REPORT_ALL(
        ch_summary_pass.concat(ch_summary_not_pass).map{meta, summaryline -> summaryline}.collect()
    )

    if(params.phylo){
        REPORT_SPECIES(
            ch_aln_stats
                .join(ch_tree, by: [0,1])
                .join(ch_dist, by: [0,1])
                .combine(REPORT_ALL.out.summary),
            file(params.microreact_template)
        )
    }

    // ---------------------------
    // Collate and save software versions
    // ---------------------------
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name    : 'CorgiSNPs_software_' + 'mqc_' + 'versions.yml',
            sort    : true,
            newLine : true
        )
        .set { ch_collated_versions }

    // ---------------------------
    // MultiQC setup
    // ---------------------------
    ch_multiqc_config = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params       = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary  = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files     = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // ---------------------------
    // MULTIQC
    // ---------------------------
    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    // Path to MultiQC HTML report
    multiqc_report = MULTIQC.out.report.toList()
    // Channel of versions.yml files from all stages
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
