process EXTRACT_REGIONS {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(ref), path(genome), path(genes)

    output:
    tuple val(meta), path("*.fastq.gz"),             emit: reads , optional:true
    tuple val(meta), path("genome_masked.fasta.gz"), emit: ref
    path "versions.yml",                             emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    tool = 'coordcutter'
    """
    ${tool} \\
        --bam ${bam} \\
        --ref ${ref} \\
        --genome ${genome} \\
        --gff ${genes} \\
        --genes "${params.amr_genes}"

    gzip genome_masked.fasta
    
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
