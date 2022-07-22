#!/usr/bin/env nextflow

// prepare lineages and representative files

workflow PREPARE_LINEAGES_REPRESENTATIVE {

    take:
    prepare_segment_indices

    main:

    if (params.segment == 'Gn') {
        fasta = 'Gn'? params.segments[ 'RVFV' ][ 'Gn' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "representative sequences file not found: ${fasta}"}
        println ch_fasta.view()


        lineages = 'Gn'? params.segments[ 'RVFV' ][ 'Gn' ].csv ?: false : false
        ch_lineages = Channel
            .fromPath(lineages, checkIfExists: true)
            .ifEmpty { exit 1, "lineages file not found: ${lineages}"}


    }

    if (params.segment == 'S') {
        fasta = 'RVFV'? params.segments[ 'RVFV' ][ 'S' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "representative sequences file not found: ${fasta}"}
        println ch_fasta.view()


        lineages = 'RVFV'? params.segments[ 'RVFV' ][ 'S' ].csv ?: false : false
        ch_lineages = Channel
            .fromPath(lineages, checkIfExists: true)
            .ifEmpty { exit 1, "lineages file not found: ${lineages}"}

    }

    if (params.segment == 'M') {
        fasta = 'RVFV'? params.segments[ 'RVFV' ][ 'M' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "representative sequences file not found: ${fasta}"}

        println ch_fasta.view()


        lineages = 'RVFV'? params.segments[ 'RVFV' ][ 'M' ].csv ?: false : false
        ch_lineages = Channel
            .fromPath(lineages, checkIfExists: true)
            .ifEmpty { exit 1, "lineages file not found: ${lineages}"}

    }

    if (params.segment == 'L') {
        fasta = 'RVFV'? params.segments[ 'RVFV' ][ 'L' ].fasta ?: false : false
        ch_fasta = Channel
            .fromPath(fasta, checkIfExists: true)
            .ifEmpty { exit 1, "representative sequences file not found: ${fasta}"}
        println ch_fasta.view()


        lineages = 'RVFV'? params.segments[ 'RVFV' ][ 'L' ].csv ?: false : false
        ch_lineages = Channel
            .fromPath(lineages, checkIfExists: true)
            .ifEmpty { exit 1, "lineages file not found: ${lineages}"}

    }

    emit:
    fasta                 = ch_fasta             // path: representative.fasta
    lineages              = ch_lineages          // path: lineages.csv
}