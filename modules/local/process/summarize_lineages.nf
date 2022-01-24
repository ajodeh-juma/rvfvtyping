// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SUMMARY_LINEAGES {
    tag "summarize_lineages"
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
    path query
    path assignment
    tuple val(meta) path(lineage)

    output:
    path "*.csv",         emit: summary

    script:
    def software  = getSoftwareName(task.process)
    // def blast     = diamond.name != 'None' ? "--blast $diamond" : ''

    """
    report.py \\
        --query $query \\
        --assignment $assignment \\
        --lineage-csv $lineage
    """
    
}
