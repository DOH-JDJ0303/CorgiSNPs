process SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/staphb/snpeff:5.2f'

    input:
    tuple val(meta), path(vcf), path(tbi), val(species)
    path snpeff_db, stageAs: 'snpEff'

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        ${species} \\
        -c "${snpeff_db}/snpEff.config" \\
        -csvStats ${prefix}.csv \\
        $args \\
        $vcf \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    touch ${prefix}.csv
    touch ${prefix}.html
    touch ${prefix}.genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

}
