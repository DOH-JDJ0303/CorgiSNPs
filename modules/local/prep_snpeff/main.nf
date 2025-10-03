process PREP_SNPEFF {
    tag "${species}"
    label 'process_low'
    
        
    input:
    val species
    path snpeff_db

    output:
    tuple val(species), path("${sp_dir}/sequences.fa"), path("${sp_dir}/genes.gff"), path("${sp_dir}/targets.json"), emit: files, optional: true

    script:
    def args = task.ext.args ?: ''
    sp_dir = "${snpeff_db}/data/${species}"
    """
    """
}
