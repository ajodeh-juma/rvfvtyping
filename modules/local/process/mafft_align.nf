// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MAFFT_ALIGN {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::mafft=7.475' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mafft:7.475--h779adbc_1'
    } else {
        container 'quay.io/biocontainers/mafft:7.475--h779adbc_1'
    }

    input:
    path fasta

    output:
    path '*.fasta'                           , emit: alignment
    path '*.log'                             , emit: log
    path '*.version.txt'                     , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix   = fasta.baseName

    // --reorder \\
    // --anysymbol \\
    //  --nomemsave \\
    //    --adjustdirection \\

    """
    mafft \\
        --thread $task.cpus \\
        $fasta 1> \\
        ${prefix}.align.fasta 2> \\
        ${prefix}.log
    
    echo \$(mafft --version 2>&1) | sed -e 's/^*.v\$//' > ${software}.version.txt
    """
}