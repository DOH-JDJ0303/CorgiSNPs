//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SHOVILL      } from '../../../modules/local/shovill/main'
include { GAMBIT_QUERY } from '../../../modules/local/gambit/main'
include { SUBTYPE      } from '../../../modules/local/subtype/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CLASSIFY {

    take:
    // ch_samplesheet: [ val(meta), path(reads) ]
    ch_samplesheet

    main:

    // Collector for version files
    ch_versions = Channel.empty()

    // -------------------------------------------------------------------------
    // MODULE: de novo assembly (produces contigs for classification)
    // -------------------------------------------------------------------------
    SHOVILL(
        ch_samplesheet
    )

    // -------------------------------------------------------------------------
    // MODULE: Species assignment via GAMBIT (uses contigs)
    // -------------------------------------------------------------------------
    GAMBIT_QUERY(
        SHOVILL.out.contigs,
        params.gambit_db,
        params.gambit_h5_dir
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions.first())

    // -------------------------------------------------------------------------
    // MODULE: Subtyping (joins assigned taxa with contigs)
    // -------------------------------------------------------------------------
    SUBTYPE(
        GAMBIT_QUERY.out.taxa.join(SHOVILL.out.contigs, by: 0),
        params.subtype_db
    )
    ch_versions = ch_versions.mix(SUBTYPE.out.versions.first())

    // -------------------------------------------------------------------------
    // Parse subtype CSV to a normalized stream:
    // [ meta, sanitized_species, sanitized_subtype ]
    // -------------------------------------------------------------------------
    SUBTYPE.out.subtype
        .splitCsv(header: true)
        .map { meta, data -> [ meta, Utils.sanitize(data['taxon']), Utils.sanitize(data['subtype']) ] }
        .set { ch_class }

    // -------------------------------------------------------------------------
    // Update the species and subtype values
    // -------------------------------------------------------------------------
    ch_samplesheet
        .join(ch_class)
        .map { meta, reads, species, subtype ->
            def new_meta = meta + [species: species ? species : 'no_species', subtype: subtype ? subtype : 'no_subtype']
            [ new_meta, reads ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    denovo      = SHOVILL.out.contigs
    species     = GAMBIT_QUERY.out.taxa
    subtype     = SUBTYPE.out.subtype
    versions    = ch_versions
}
