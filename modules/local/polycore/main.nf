process POLYCORE {
    tag "${prefix}"
    label 'process_high'
        
    input:
    tuple val(species), val(subtype), path(assemblies), path(ref)

    output:
    tuple val(species), val(subtype), path("core.aln"),      emit: snps
    tuple val(species), val(subtype), path("core.full.aln"), emit: full
    tuple val(species), val(subtype), path("summary.csv"),   emit: csv
    tuple val(species), val(subtype), path("fconst.txt"),    emit: fconst
    tuple val(species), val(subtype), path("dist_long.csv"), emit: dist_long
    tuple val(species), val(subtype), path("dist_wide.csv"), emit: dist_wide
    tuple val(species), val(subtype), path("*.html"),        emit: plot, optional: true
    path "versions.yml",                                     emit: versions


    script:
    def args = task.ext.args ?: ''
    prefix = "${species}-${subtype}"
    """
    polycore \\
        --ref ${ref} \\
        --sample ${assemblies} \\
        --min-gf ${params.min_genome_fraction} \\
        --min-cf ${params.min_core_fraction} \\
        --ploidy ${params.ploidy} \\
        ${args}
        
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polycore: "\$(polycore --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
