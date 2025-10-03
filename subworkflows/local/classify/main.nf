import Utils

//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SHOVILL      } from '../../../modules/nf-core/shovill/main'
include { GAMBIT_QUERY } from '../../../modules/local/gambit/main'
include { SUBTYPE      } from '../../../modules/local/subtype/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CLASSIFY {

    take:
    // ch_reads: [ val(meta), path(reads) ]
    ch_reads
    // ch_meta : [ val(meta), val(species), val(subtype), path(reference) ]
    ch_meta

    main:

    // Collector for version files
    ch_versions = Channel.empty()

    // -------------------------------------------------------------------------
    // MODULE: de novo assembly (produces contigs for classification)
    // -------------------------------------------------------------------------
    SHOVILL(
        ch_reads
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
    // Prefer provided species/subtype; otherwise fall back to derived values. (Shouldn't be needed but still included)
    // Preserve original reference path.
    // Output shape: [ meta, species, subtype, reference ]
    // -------------------------------------------------------------------------
    ch_meta
        .join(ch_class)
        .map { meta, species, subtype, reference, derived_species, derived_subtype ->
            [ meta,
              species  ? species  : derived_species,
              subtype  ? subtype  : derived_subtype,
              reference ]
        }
        .set { ch_meta }

    // Final catch-all
    ch_meta.map{ meta, sp, sb, ref -> [meta, sp ? sp : 'no_species', sb ? sb : 'no_subtype', ref ] }.set{ch_meta}

    emit:
    meta     = ch_meta
    denovo   = SHOVILL.out.contigs
    species  = GAMBIT_QUERY.out.taxa
    subtype  = SUBTYPE.out.subtype
    versions = ch_versions
}
