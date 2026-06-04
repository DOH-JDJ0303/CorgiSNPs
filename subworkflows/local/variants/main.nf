//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN_READS        } from '../../../modules/local/align_reads/main'
include { SAMTOOLS_MPILEUP   } from '../../../modules/nf-core/samtools/mpileup/main'
include { LOWSITES           } from '../../../modules/local/lowsites/main'
include { FREEBAYES          } from '../../../modules/local/freebayes/main'
include { FILTER_VCF         } from '../../../modules/local/filter_vcf/main'
include { CREATE_MASK        } from '../../../modules/local/create_mask/main'
include { BCFTOOLS_CONSENSUS } from '../../../modules/local/bcftools/consensus/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VARIANTS {

    take:
    // ch_reads: [ val(meta), path(reads) ]
    ch_samplesheet
    // make_consensus: boolean flag controlling consensus generation
    make_consensus

    main:

    // Collectors
    ch_versions = Channel.empty()
    ch_aln      = Channel.empty()

    // -------------------------------------------------------------------------
    // Load reference DB JSON and fan-out entries
    // -------------------------------------------------------------------------
    Channel
        .fromPath(params.reference_db)
        .splitJson()
        .map { rec ->
            if (!(rec instanceof Map)) {
                throw new IllegalArgumentException("Each JSON item must be an object/map. Got: ${rec?.getClass()?.name}")
            }

            // Required keys present & non-empty?
            def required = ['species','subtype','reference','ploidy']
            def missing = required.findAll { k ->
                !rec.containsKey(k) || rec[k] == null || (
                    rec[k] instanceof Collection ? rec[k].isEmpty() : (rec[k] instanceof String && rec[k].trim().isEmpty())
                )
            }
            if (missing) {
                throw new IllegalArgumentException(
                    "Invalid record in '${params.reference_db}': missing/empty ${missing.join(', ')}.\nRecord: ${groovy.json.JsonOutput.toJson(rec)}"
                )
            }

            // Normalize fields (allow string or list for species/subtype)
            def species  = (rec.species  instanceof Collection) ? rec.species  : [rec.species]
            def subtype  = (rec.subtype  instanceof Collection) ? rec.subtype  : [rec.subtype]
            def refPath  = rec.reference as String
            def ploidy   = rec.ploidy

            // Enforce .gz and file existence
            if (!refPath.endsWith('.gz')) {
                throw new IllegalArgumentException("Reference must be a .gz file. Got: '${refPath}'")
            }
            def refFile = file("${projectDir}/assets/").resolve(refPath)
            if (!refFile.exists()) {
                throw new IllegalArgumentException("Reference file not found: '${refFile}'")
            }

            // Return a clean, predictable shape
            [ species: species, subtype: subtype, reference: refFile, ploidy: ploidy ]
        }
        .set { ch_refs_db }

    // -------------------------------------------------------------------------
    // For 'auto' rows: find matching reference by species/subtype; then unify with 'manual'
    // Output ch_meta retains shape: [ meta, sp, sb, ref ]
    // -------------------------------------------------------------------------
    ch_samplesheet = ch_samplesheet
        .filter{ meta, reads -> ! meta.reference }
        .combine(ch_refs_db)
        .map { meta, reads, refData ->
            def supplied     = meta.reference && meta.ploidy ? true : false
            def speciesMatch = refData.species.any  { Utils.sanitize(it) == Utils.sanitize(meta.species) }
            def subtypeMatch = refData.subtype.any { Utils.sanitize(it) == Utils.sanitize(meta.subtype) }
            def match        = supplied ? false : speciesMatch && subtypeMatch
            def reference    = match ? refData.reference : meta.reference 
            def ploidy       = match ? refData.ploidy : meta.ploidy
            def new_meta     = meta + [reference: reference, ploidy: ploidy]

            return [new_meta, reads]
        }
        .filter{ meta, reads -> meta.reference }
        .concat(
            ch_samplesheet.filter{ meta, reads -> meta.reference }
        )

    // -------------------------------------------------------------------------
    // ALIGN_READS
    // -------------------------------------------------------------------------
    ALIGN_READS(
        ch_samplesheet.map{ meta, reads -> [meta, reads, meta.reference] }
    )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions.first())
    ALIGN_READS.out.mapped.set { ch_mapped }

    // -------------------------------------------------------------------------
    // FREEBAYES (variant calling)
    // -------------------------------------------------------------------------
    FREEBAYES(
        ch_mapped.map{ meta, bam, bai -> [meta, bam, bai, meta.reference] }
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    // -------------------------------------------------------------------------
    // FILTER_VCF (post-calling filtering)
    // -------------------------------------------------------------------------
    FILTER_VCF(
        FREEBAYES.out.vcf
    )
    ch_versions = ch_versions.mix(FILTER_VCF.out.versions.first())

    // -------------------------------------------------------------------------
    // Optional: consensus generation + (optional) push to DB
    // When enabled, compute mpileup/low-sites for masking and build consensus.
    // -------------------------------------------------------------------------
    ch_depth = ch_samplesheet.map{[it[0], []]}
    if (make_consensus) {

        // SAMTOOLS mpileup (depth & pileup summaries)
        // mpileup inputs: [ meta, bam, opts ]; opts left as [] to preserve behavior
        SAMTOOLS_MPILEUP(
            ch_mapped.map { meta, bam, bai -> [ meta, bam, [] ] },
            [ null, [] ]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

        // LOWSITES (derive low-coverage sites and summary)
        LOWSITES(
            SAMTOOLS_MPILEUP.out.mpileup
        )
        ch_versions = ch_versions.mix(LOWSITES.out.versions.first())
        ch_depth    = LOWSITES.out.summary

        CREATE_MASK (
            LOWSITES.out.bed
                .join(FILTER_VCF.out.bed)
        )

        // Build consensus using filtered SNVs + reference + low-site mask
        BCFTOOLS_CONSENSUS(
            FILTER_VCF
                .out
                .snvs
                .map{ meta, vcf, tbi -> [meta, vcf, tbi, meta.reference] }
                .join(CREATE_MASK.out.bed)

        )
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())
        ch_aln      = BCFTOOLS_CONSENSUS.out.fasta

        // Optional: copy consensus into a DB layout by species/subtype
        if (params.push) {
            BCFTOOLS_CONSENSUS.out.fasta
                .subscribe { meta, fa ->
                    fa.copyTo(
                        file(params.db)
                            .resolve(meta.species)
                            .resolve(meta.subtype)
                            .resolve(fa.name)
                    )
                }
        }
    }

    emit:
    samplesheet = ch_samplesheet
    bam         = ch_mapped
    vcf         = FILTER_VCF.out.filt
    aln         = ch_aln
    depth       = ch_depth
    versions    = ch_versions
}
