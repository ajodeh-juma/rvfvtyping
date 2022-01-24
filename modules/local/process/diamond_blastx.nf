// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::diamond=2.0.9' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container 'quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0'
    }

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.blastx.txt')  , emit: txt
    path '*.version.txt'                   , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/.dmnd//'`

    diamond blastx \\
        --threads $task.cpus \\
        --query $fasta \\
        --db \$DB \\
        --out ${prefix}.dmnd.blastx.txt \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle
    
    echo \$(diamond version 2>&1) | tr -d [a-z] | sed -e 's/^[[:space:]]*//' > ${software}.version.txt
    """
}