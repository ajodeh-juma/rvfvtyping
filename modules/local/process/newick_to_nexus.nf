// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BIOPYTHON_NEWICK_TO_NEXUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::biopython=1.78" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/biopython:1.78"
    } else {
        container "quay.io/biocontainers/biopython:1.78"
    }

    input:
    tuple val(meta), path(in_tree)

    output:
    tuple val(meta), path("*.nexus")         , emit: out_tree
    path '*.version.txt'                     , emit: version

    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def in_format   = "newick"
    def out_format  = "nexus"
    
    """
    convertTreeFormats.py \\
        --in-tree $in_tree \\
        --in-format ${in_format} \\
        --out-tree ${prefix}.nexus \\
        --out-format ${out_format}

    echo \$(python -c 'import Bio; print(Bio.__version__)') > ${software}.version.txt
    """
}
