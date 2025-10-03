import Utils

//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { POLYCORE } from '../../../modules/local/polycore/main'
include { IQTREE   } from '../../../modules/local/iqtree/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PHYLO {

    take:
    // ch_aln : [ val(meta), path(aln) ]  // consensus alignments (label typo kept)
    ch_aln
    // ch_meta: [ val(meta), val(species), val(subtype), path(reference) ]
    ch_meta

    main:

    // Collector for version files
    ch_versions = Channel.empty()

    // -------------------------------------------------------------------------
    // Join alignments with meta, group by (species, subtype), and assemble a
    // sample set consisting of current alignments plus any cached consensus
    // sequences from params.db not already present in the current set.
    // -------------------------------------------------------------------------
    ch_aln
        .join(ch_meta, by: 0)
        .groupTuple(by: [2, 3]) // group by species (idx 2) and subtype (idx 3)
        .map { meta, aln, sp, sb, refs ->
            // Current sample names (filenames) from the alignment paths
            def new_aln = aln.collect { it.getName() }

            // Ensure a single unique reference for (sp, sb)
            def ref = refs.unique()
            if (ref.size() > 1) {
                exit("ERROR: Multiple references associated with species: ${sp} subtype: ${sb}\n${ref}")
            }

            // Collect cached consensus sequences for this (sp/sb) from params.db
            def sp_sb_db       = file(params.db).resolve(sp).resolve(sb)
            def sp_db_cache    = sp_sb_db.exists() ? sp_sb_db.list().collect { sp_sb_db.resolve(it) } : []
            def filtered_cache = sp_db_cache.findAll { !new_aln.contains(it.getName()) }

            // Combine current alignments with cache, keeping ref as a single path
            def samples = aln + filtered_cache
            [ sp, sb, samples, ref[0] ]
        }
        .set { ch_sp_sb }
    ch_sp_sb

    // -------------------------------------------------------------------------
    // POLYCORE: derive SNP alignment, distance matrix, and per-sample stats
    // Input: [ species, subtype, samples, reference ]
    // -------------------------------------------------------------------------
    POLYCORE(
        ch_sp_sb
    )
    ch_versions = ch_versions.mix(POLYCORE.out.versions.first())

    // Summarize per-(species, subtype) sample counts for downstream filtering
    POLYCORE.out.csv
        .splitCsv(header: true)
        .groupTuple(by: [0, 1])      // group by species, subtype
        .map { sp, sb, data -> [ sp, sb, data.size() ] }
        .set { ch_sample_counts }

    // -------------------------------------------------------------------------
    // IQTREE: build trees only when there are >2 samples for the clade
    // Joins: SNP alignment + constant-site file + sample count
    // -------------------------------------------------------------------------
    IQTREE(
        POLYCORE.out.snps
            .join(POLYCORE.out.fconst, by: [0, 1])
            .join(ch_sample_counts,    by: [0, 1])
            .filter { sp, sb, aln, const_sites, count -> count > 2 }
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions.first())

    emit:
    meta      = ch_meta
    aln_stats = POLYCORE.out.csv
    tree      = IQTREE.out.tree
    dist      = POLYCORE.out.dist_wide
    versions  = ch_versions
}
