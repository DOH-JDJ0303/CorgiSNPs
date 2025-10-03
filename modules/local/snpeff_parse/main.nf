process SNPEFF_PARSE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(targets_json)

    output:
    tuple val(meta), path("variants_full.csv"),    emit: full, optional: true
    tuple val(meta), path("variants_targets.csv"), emit: target, optional: true
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    tool = 'snpeff_parser.py'
    """
    ${tool} \\
        -r ${targets_json} \\
        ${vcf}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """

}
