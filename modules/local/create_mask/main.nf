process CREATE_MASK {
    tag "${meta.id}"
    label 'process_low'
        
    input:
    tuple val(meta), path(lowsites_bed), path(vcf_beds)

    output:
    tuple val(meta), path("*.final-mask.bed"), emit: bed
    path "versions.yml",                       emit: versions
    

    script:
    def args = task.ext.args ?: '' 
    """
    # Sort + merge
    cat ${lowsites_bed} *.filt_vcf.bed | cut -f 1-3 | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge > rough_mask.bed
    cat *.pass_snp.bed | LC_ALL=C sort -k1,1 -k2,2n -k3,3n | bedtools merge > keep.bed

    # Subtract PASS SNP positions from the mask so we don't mask true passing SNPs
    bedtools subtract -a rough_mask.bed -b keep.bed > ${meta.id}.final-mask.bed

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sha256sum: "\$(bedtools --version | sed -n '1p' | cut -f 2 -d ' ')"
    END_VERSIONS
    """
}
