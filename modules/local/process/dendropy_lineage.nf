// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DENDROPY_LINEAGE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::dendropy=4.5.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/dendropy:4.5.2--pyh3252c3a_0"
    } else {
        container "quay.io/biocontainers/dendropy:4.5.2--pyh3252c3a_0"
    }

    input:
    tuple val(meta), path(in_tree)

    output:
    tuple val(meta), path("*.csv")           , emit: lineage_csv
    path '*.version.txt'                     , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def collapse    = 0.000005

    """
    
    assignLineage.py  \\
        --separator '|' \\
        --index 1 \\
        --taxon ${prefix} \\
        --collapse_to_polytomies ${collapse} \\
        --input $in_tree \\
        --output ${prefix}.lineage.csv

    
    echo \$(python -c 'import dendropy; print(dendropy.__version__)') > ${software}.version.txt
    """
}
