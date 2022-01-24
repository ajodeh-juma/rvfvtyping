// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GUNZIP_DATABASE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? 'conda-forge::curl=7.76.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/curl:7.45.0--1'
    } else {
        container 'quay.io/biocontainers/curl:7.45.0--1'
    }

    input:
    tuple val(meta), path(db)

    output:
    path 'diamond'                 , emit: database
    path '*.version.txt'           , emit: version

    script:
    def software        = getSoftwareName(task.process)
    def lastPath        = "${db}".lastIndexOf(File.separator)
    def lastExt         = "${db}".lastIndexOf(".")
    def base            = "${db}".substring(lastPath+1,lastExt)
    
    """
    mkdir diamond
    gunzip --keep -f -c $db > diamond/$base
    echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//' > ${software}.version.txt
    """
}