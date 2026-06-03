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
    // ch_samplesheet: [ meta, reads ]
    ch_samplesheet
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
        ch_samplesheet.map { meta, reads -> meta.species }.unique(),
        params.snpeff_db
    )

    // Tie SnpEff assets (fa/gff/json) back to each sample's meta/ref
    ch_samplesheet
        .map    { meta, reads -> [ meta.species, meta ] }
        .combine(PREP_SNPEFF.out.files, by: 0)
        .map    { sp, meta, fa, gff, json -> [ meta, fa, gff, json ] }
        .set    { ch_snpeff_files }

    // ---------------------------------------------------------------------
    // Compare sample reference vs SnpEff FASTA (checksum-based)
    // ---------------------------------------------------------------------
    COMPARE_REFS(
        ch_snpeff_files.map { meta, fa, gff, json -> [ meta, meta.reference, fa ] }
    )
    ch_versions = ch_versions.mix(COMPARE_REFS.out.versions.first())

    // Determine per-sample whether all checksums agree (same reference)
    ch_ref_match = COMPARE_REFS.out.txt
        .splitText()
        .map       { meta, line -> [ meta, line.split()[0] ] } // [meta, checksum]
        .groupTuple()
        .map       { meta, checksums -> [ meta, checksums.unique().size() == 1 ] }

    // Branch: samples whose references match SnpEff vs those that don't
    ch_ref_match
        .branch { meta, matches ->
            same_ref:  matches
            diff_ref: !matches
        }
        .set { ch_ref_branches }

    // ---------------------------------------------------------------------
    // For differing references: extract target regions and re-call variants
    // ---------------------------------------------------------------------
    EXTRACT_REGIONS(
        ch_ref_branches.diff_ref
            .map  { meta, matches -> meta } // strip boolean
            .join(
                ch_bam.map { meta, bam, bai -> [ meta, bam ] }
            )
            .join(
                ch_snpeff_files.map { meta, fa, gff, json -> [ meta, meta.reference, fa, gff ] }
            )
    )
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions.first())
    ch_regions  = EXTRACT_REGIONS
        .out
        .results
        .map{ meta, reads, ref ->
                def new_meta = meta + [reference: ref]
                return [new_meta, reads]
        }
    // Call variants on extracted reads against the extracted reference
    VARIANTS(
        ch_regions,
        false
    )
    ch_versions = ch_versions.mix(VARIANTS.out.versions.first())
    ch_new_vcf  = VARIANTS.out.vcf

    // ---------------------------------------------------------------------
    // SnpEff: annotate either the original VCF (same_ref) or re-called VCF
    // ---------------------------------------------------------------------
    ch_for_snpeff = ch_ref_branches
        .same_ref
        .map   { meta, matches -> meta }
        .join  (ch_vcf)
        .concat(ch_new_vcf)

    SNPEFF(
        ch_for_snpeff,
        params.snpeff_db
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())

    // Parse SnpEff results using species JSON emitted by PREP_SNPEFF
    // Match annotated VCFs to species JSON by stable sample ID (meta.id), not full meta
    SNPEFF_PARSE(
        SNPEFF
            .out
            .vcf
            .map { meta, vcf -> [ meta.id, meta, vcf ] }
            .join(
                ch_snpeff_files
                    .map { meta, fa, gff, json -> [ meta.id, json ] },
                by: 0
            )
            .map { id, meta, vcf, json -> [ meta, vcf, json ] }
    )
    ch_versions = ch_versions.mix(SNPEFF_PARSE.out.versions.first())

    emit:
    summary  = SNPEFF_PARSE.out.target
    versions = ch_versions
}
