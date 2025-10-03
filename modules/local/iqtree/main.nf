process IQTREE {
    tag "${prefix}"
    label 'process_high'
    
    input:
    tuple val(species), val(subtype), path(aln), path(const_sites), val(count)

    output:
    tuple val(species), val(subtype), path("*.nwk"), emit: tree
    path 'versions.yml',                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args         = task.ext.args ?: ''
    prefix       = "${species}-${subtype}"
    bootstrap    = count > 4 ? '-B 1000' : ''
    tree_ext     = count > 4 ? 'contree' : 'treefile'
    """
    # run IQTREE2
    iqtree2 \\
        -s ${aln} \\
        -fconst \$(cat ${const_sites}) \\
        -T ${task.cpus} \\
        ${args} \\
        ${bootstrap}

    mv *.${tree_ext} ${prefix}.nwk

    #### VERSION INFO ####
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(iqtree2 --version | head -n 1 | cut -f 4  -d ' ')
    END_VERSIONS
    """
}