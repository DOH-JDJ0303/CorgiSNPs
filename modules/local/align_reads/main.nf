process ALIGN_READS {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9ac054213e67b3c9308e409b459080bbe438f8fd6c646c351bc42887f35a42e7/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd' }"

    input:
    tuple val(meta), path(reads), path(ref)

    output:
    tuple val(meta), path("*.mapped.bam"), path("*.mapped.bam.bai"),     emit: mapped , optional:true
    tuple val(meta), path("*.unmapped.bam"), path("*.unmapped.bam.bai"), emit: unmapped , optional:true
    path  "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bwa-mem2 index ${ref}
    bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        ${ref} \\
        $reads \\
        > ${meta.id}.all.bam

    # Extract mapped reads
    samtools view -b -F 4 ${meta.id}.all.bam | samtools sort > ${meta.id}.mapped.bam
    samtools index ${meta.id}.mapped.bam
    
    # Extract unmapped reads
    samtools view -b -f 4 ${meta.id}.all.bam | samtools sort > ${meta.id}.unmapped.bam
    samtools index ${meta.id}.unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
