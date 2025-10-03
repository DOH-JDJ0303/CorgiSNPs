process COMPARE_REFS {
    tag "${meta.id}"
    label 'process_low'
        
    input:
    tuple val(meta), path(ref1), path(ref2)

    output:
    tuple val(meta), path("sha256sum.txt"), emit: txt
    path "versions.yml",                    emit: versions
    

    script:
    def args = task.ext.args ?: '' 
    """
    sha256sum ${ref1} ${ref2} > sha256sum.txt

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sha256sum: "\$(sha256sum --version | sed -n '1p' | cut -f 4 -d ' ')"
    END_VERSIONS
    """
}
