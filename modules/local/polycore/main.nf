process POLYCORE {
    tag "${prefix}"
    label 'process_high'
        
    input:
    tuple val(meta), path(ref), path(assemblies)

    output:
    tuple val(meta), path("core.aln"),      emit: snps
    tuple val(meta), path("core.full.aln"), emit: full
    tuple val(meta), path("summary.csv"),   emit: csv
    tuple val(meta), path("full.csv"),      emit: full_csv
    tuple val(meta), path("fconst.txt"),    emit: fconst
    tuple val(meta), path("dist_long.csv"), emit: dist_long
    tuple val(meta), path("dist_wide.csv"), emit: dist_wide
    tuple val(meta), path("*.html"),        emit: plot, optional: true
    path "versions.yml",                    emit: versions


    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.species}-${meta.subtype}"
    """
    polycore \\
        ${ref} ${assemblies} \\
        --min-gf ${params.min_genome_fraction} \\
        --min-cf ${params.min_core_fraction} \\
        --ploidy ${meta.ploidy} \\
        ${args}
        
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polycore: "\$(polycore --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
