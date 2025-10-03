process SUMMARYLINE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(read_stats), path(denovo), path(species), path(subtype)
    path(samplesheet)
    path(ncbi_stats)

    output:
    tuple val(meta), path("*-summary.csv"), emit: summary
    path "versions.yml",                    emit: versions

    script:
    def tool = 'summaryline.py'
    def args = []
    if (read_stats) args << "--read_stats \"${read_stats}\""
    if (denovo)     args << "--denovo \"${denovo}\""
    if (species)    args << "--species \"${species}\""
    if (subtype)    args << "--subtype \"${subtype}\""

    """
    ${tool} \\
        --sample "${meta.id}" \\
        --samplesheet ${samplesheet} \\
        --ncbi_stats ${ncbi_stats} \\
        --min_depth ${params.min_depth_qc} \\
        --min_qual ${params.min_q30_rate_qc} \\
        --max_z_score ${params.max_z_score_qc} \\
        ${args.join(' ')}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}

process ADD_AMR {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(summaryline), path(amr)

    output:
    tuple val(meta), path("*-summary.csv", includeInputs: true), emit: summary
    path "versions.yml",                    emit: versions

    script:
    def tool = 'add_amr.py'
    def args = []

    """
    ${tool} \\
        --sample "${meta.id}" \\
        --summaryline "${summaryline}" \\
        --amr ${amr} \\
        ${args.join(' ')}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}

process REPORT_ALL {
    label 'process_low'
        
    input:
    path summarylines

    output:
    path "MycoSNP-summary.csv", emit: summary
    path "versions.yml",        emit: versions
    
    script:
    tool = 'report_all.py'
    """
    ${tool} ${summarylines}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}

process REPORT_SPECIES {
    tag "${prefix}"
    label 'process_low'
        
    input:
    tuple val(species), val(subtype), path(aln_stats), path(tree), path(dist), path(summary)
    path microreact_template

    output:
    tuple val(species), val(subtype), path("*"), emit: results
    path "versions.yml",                         emit: versions
    
    script:
    prefix = "${species}-${subtype}"
    tool = 'report_species.py'
    """
    ${tool} \\
        --prefix ${prefix} \\
        --aln_stats "${aln_stats}" \\
        --dist "${dist}" \\
        --summary ${summary} \\
        --tree "${tree}" \\
        --strong_link ${params.strong_link_threshold} \\
        --inter_link ${params.inter_link_threshold} \\
        --partition_distance ${params.partition_distance} \\
        --microreact "${microreact_template}"

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}

