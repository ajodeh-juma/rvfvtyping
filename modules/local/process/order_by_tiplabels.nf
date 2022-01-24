// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ORDER_BY_TIPLABELS {
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
    tuple val(meta), path(tree)
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.fasta")       , emit: fasta

    script:
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    orderseqs_by_tiplabels.py \\
        --alignment $alignment \\
        --tree $tree \\
        --outfile ${prefix}.orderedseq.fasta
    """
}
