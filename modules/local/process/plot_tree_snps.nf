// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLOT_TREE_SNPS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-argparse=2.0.3 conda-forge::r-dplyr=1.0.5 conda-forge::r-ggplot2=3.3.4 conda-forge::r-rcolorbrewer=1.1_2 conda-forge::r-rvcheck=0.1.8 conda-forge::r-tidytree=0.3.4 bioconda::bioconductor-biostrings=2.58.0 bioconda::bioconductor-ggtree=2.4.1 bioconda::bioconductor-treeio=1.14.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-base:3.5.1"
    } else {
        container "quay.io/biocontainers/r-base:3.5.1"
    }

    input:
    tuple val(meta), path(tree)
    tuple val(meta), path(snps)
    path (lineages)

    output:
    tuple val(meta), path("*.pdf")      , emit: pdf

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    plot_tree_snps.r \\
        --treefile $tree \\
        --snps $snps \\
        --metadata $lineages \\
        --prefix ${prefix}
    """
}