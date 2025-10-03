import Utils

//
// Subworkflow for AMR analysis in the mycosnp pipeline
//

include { PREP_SNPEFF     } from '../../../modules/local/prep_snpeff/main'
include { COMPARE_REFS    } from '../../../modules/local/compare_refs/main'
include { EXTRACT_REGIONS } from '../../../modules/local/extract_regions/main'
include { VARIANTS        } from '../../../subworkflows/local/variants'
include { SNPEFF          } from '../../../modules/local/snpeff/main'
include { SNPEFF_PARSE    } from '../../../modules/local/snpeff_parse/main'


workflow AMR {

    take:
    // ch_meta: [ meta, species, subtype, reference ]
    ch_meta
    // ch_bam : [ meta, bam, bai ]
    ch_bam
    // ch_vcf : [ meta, vcf, csi ]
    ch_vcf

    main:
    // Collector for version files
    ch_versions = Channel.empty()

    // ---------------------------------------------------------------------
    // Prepare SnpEff reference files once per species in the batch
    // ---------------------------------------------------------------------
    PREP_SNPEFF(
        ch_meta.map { meta, sp, sb, ref -> sp }.unique(),
        params.snpeff_db
    )

    // Tie SnpEff assets (fa/gff/json) back to each sample's meta/ref
    ch_meta
        .map    { meta, sp, sb, ref -> [ sp, meta, ref ] }
        .combine(PREP_SNPEFF.out.files, by: 0)
        .map    { sp, meta, ref, fa, gff, json -> [ meta, sp, ref, fa, gff, json ] }
        .set    { ch_snpeff_files }

    // ---------------------------------------------------------------------
    // Compare sample reference vs SnpEff FASTA (checksum-based)
    // ---------------------------------------------------------------------
    COMPARE_REFS(
        ch_snpeff_files.map { meta, sp, ref, fa, gff, json -> [ meta, ref, fa ] }
    )
    ch_versions = ch_versions.mix(COMPARE_REFS.out.versions.first())

    // Determine per-sample whether all checksums agree (same reference)
    ch_ref_match = COMPARE_REFS.out.txt
        .splitText()
        .map       { meta, line -> [ meta, line.split()[1] ] } // [meta, checksum]
        .groupTuple()
        .map       { meta, checksums -> [ meta, checksums.unique().size() == 1 ] }

    // Branch: samples whose references match SnpEff vs those that don't
    ch_ref_match
        .branch { meta, matches ->
            same_ref:  matches
            diff_ref: !matches
        }
        .set { ref_branches }

    // ---------------------------------------------------------------------
    // For differing references: extract target regions and re-call variants
    // ---------------------------------------------------------------------
    EXTRACT_REGIONS(
        ref_branches.diff_ref
            .map  { meta, matches -> meta } // strip boolean
            .join(ch_bam.map { meta, bam, bai -> [ meta, bam ] })
            .join(ch_snpeff_files.map { meta, sp, ref, fa, gff, json -> [ meta, ref, fa, gff ] })
    )
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions.first())

    // Call variants on extracted reads against the extracted reference
    VARIANTS(
        EXTRACT_REGIONS.out.reads,
        EXTRACT_REGIONS.out.ref.map { meta, ref -> [ meta, [], [], ref ] },
        false
    )
    ch_versions = ch_versions.mix(VARIANTS.out.versions.first())
    VARIANTS.out.vcf.set { ch_new_vcf }

    // ---------------------------------------------------------------------
    // SnpEff: annotate either the original VCF (same_ref) or re-called VCF
    // ---------------------------------------------------------------------
    ref_branches.same_ref
        .map   { meta, matches -> meta }
        .join  (ch_vcf)
        .concat(ch_new_vcf)
        .join  (ch_meta.map { meta, species, subtype, ref -> [ meta, species ] })
        .set   { ch_for_snpeff }

    SNPEFF(
        ch_for_snpeff,
        params.snpeff_db
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())

    // Parse SnpEff results using species JSON emitted by PREP_SNPEFF
    SNPEFF_PARSE(
        SNPEFF.out.vcf.join( ch_snpeff_files.map { meta, sp, ref, fa, gff, json -> [ meta, json ] } )
    )
    ch_versions = ch_versions.mix(SNPEFF_PARSE.out.versions.first())

    emit:
    summary  = SNPEFF_PARSE.out.target
    versions = ch_versions
}
