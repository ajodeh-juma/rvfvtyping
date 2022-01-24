// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FIND_SNPS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.6.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.6.1"
    } else {
        container "quay.io/biocontainers/python:3.6.1"
    }

    input:
    tuple val(meta), path(alignment)
    tuple val(meta), path(lineages)

    output:
    tuple val(meta), path("*_seqs.csv")    , emit: representative
    tuple val(meta), path("*_snps.csv")    , emit: defining
    tuple val(meta), path("*_mask.csv")    , emit: mask

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    findSNPs.py \\
        --alignment $alignment \\
        --lineage-csv $lineages \\
        --representative-seqs-out ${prefix}.representative_seqs.csv \\
        --defining-snps-out ${prefix}.defining_snps.csv \\
        --mask-out ${prefix}.to_mask.csv

    """
}
