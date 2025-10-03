process SUBTYPE {
    tag "$meta.id"
    label 'process_low'
        
    input:
    tuple val(meta), path(gambit_results), path(seq)
    path subtype_db

    output:
    tuple val(meta), path("*_subtype.csv"), emit: subtype
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    tool = 'subtyper.py'
    """
    ${tool} \\
         --sample ${prefix} \\
         --gambit ${gambit_results} \\
         --db ${subtype_db} \\
         --seq ${seq}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
