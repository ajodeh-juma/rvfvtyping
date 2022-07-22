////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

params.summary_params = [:]

// Validate input parameters

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// check given segment
def segmentsList = ['Gn', 'S', 'M', 'L']
if (!segmentsList.contains(params.segment)) {
    exit 1, "Invalid segment option: ${params.segment}. Valid options: ${segmentsList.join(', ')}"
} else {
    ch_prefix = params.segment + '-Segment'
}

// Stage required files
if (params.segment == 'Gn') {
    fasta = file("$projectDir/segments/Gn/partial-M.representative.fasta", checkIfExists: true)
    //fasta = file("$projectDir/segments/Gn/Gn.align.trim.fasta", checkIfExists: true)
    lineages = file("$projectDir/segments/Gn/partial-M.Lineages.csv", checkIfExists: true)
    n_threshold = 90

}

if (params.segment == 'S') {
    fasta = file("$projectDir/segments/S/complete-S.align.trim.fasta", checkIfExists: true)
    lineages = file("$projectDir/segments/S/complete-S.Lineages.csv", checkIfExists: true)
    n_threshold = 90
}

if (params.segment == 'M') {
    fasta = file("$projectDir/segments/M/complete-M.align.trim.fasta", checkIfExists: true)
    lineages = file("$projectDir/segments/M/complete-M.Lineages.csv", checkIfExists: true)
    n_threshold = 90
}

if (params.segment == 'L') {
    fasta = file("$projectDir/segments/L/complete-L.align.trim.fasta", checkIfExists: true)
    lineages = file("$projectDir/segments/L/complete-L.Lineages.csv", checkIfExists: true)
    n_threshold = 90
}


if (params.input) {
    ch_fasta_files = Channel
        .fromPath(params.input)
        .map { row -> 
            def meta = [:]
            fasta = file(row, checkIfExists: true)
            records = fasta.splitFasta(record: [id: true])
            meta.id = records[0].id
            [ meta ,  fasta  ]
            }
        .ifEmpty { exit 1, 'input fasta file(s) not specified!' }
}




// diamond database file
ch_database = Channel.fromPath(params.db)
    .map { fn -> 
           def meta = [:]
           meta.id = fn.baseName
           path = fn.parent
           db = [ meta, fn] 
           return  db 
         }



// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()


// functions 
// read fasta files to get files base names
def readIn(filename) {
    Channel.from(filename)
    .map { row ->
        def meta = [:]
        def id = filename.baseName
        meta.id = "$id"
        file_info = [ meta, [filename] ]
        return file_info
    }    
}

// get percentage of Ns from summary file
def get_perc_ns_from_txt(txt) {
    def perc_ns = 0
    txt.eachLine { line ->
        if (line.contains("perc_Ns")) {
            perc_ns = line.tokenize().get(1).toFloat()
        }
    }
    return perc_ns
}

// function to filter on percent Ns
c_reset = params.monochrome_logs ? '' : "\033[0m";
c_green = params.monochrome_logs ? '' : "\033[0;32m";
c_red   = params.monochrome_logs ? '' : "\033[0;31m";


def check_input(sample, txt, ns=50.0) {
    def pass_input = [:]
    def fail_input = [:]
    perc_ns = get_perc_ns_from_txt(txt)
    if (perc_ns > ns) {
        log.info ">${c_red}>>>> $sample FAILED PERCENTAGE Ns THRESHOLD: ${perc_ns} > ${ns}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_input[sample] = perc_ns
        return false
    } else {
        pass_input[sample] = perc_ns
        return true
    }
}

// get subject id from classification output summary file
def get_sseqid(txt) {
    def sseqid = null
    txt.eachLine { line ->
        if (line.contains("sseqid")) {
            sseqid = line.tokenize().get(1).toString()
        }
    }
    return sseqid
}

// function to filter on viral segment
def check_sseqid(sample, txt, n='YP_003848705.1') {
    def pass_q = [:]
    def fail_q = [:]
    sseqid = get_sseqid(txt)
    if (sseqid != n) {
        log.info ">${c_red}>>>> $sample SAMPLE CLASSIFICATION: ${sseqid} != ${n}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_q[sample] = sseqid
        return false
    } else {
        pass_q[sample] = sseqid
        return true
    }
}

// get assignable 
def get_assignable(txt) {
    def boolean n = false
    def boolean to_keep = false
    txt.eachLine { line ->
        if (line.contains("to_keep")) {
            to_keep = line.tokenize().get(2).toBoolean()
            return to_keep
        } else {
            return to_keep
        }
    }
}


def get_tokeep(txt) {
    def int n = 0
    def int to_keep = 0
    txt.eachLine { line ->
        if (line.contains("to_keep")) {
            to_keep = line.tokenize().get(1).toInteger()
            return to_keep
        } else {
            return to_keep
        }
    }
}


////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_SOFTWARE_VERSIONS                                  } from '../modules/local/process/get_software_versions'             addParams( options: [publish_files : ['csv':'']]                 )
include { PREPARE_LINEAGES_REPRESENTATIVE                        } from '../modules/local/subworkflow/prepare_lineages_reps'         addParams( options: [:]                                          )
include { FILTER_INPUT                                           } from '../modules/local/process/filter_input'                      addParams( options: modules['filter']                            )
include { GUNZIP_DATABASE                                        } from '../modules/local/process/gunzip_database'                   addParams( options: [:]                                          )
include { DIAMOND_BLASTX                                         } from '../modules/local/process/diamond_blastx'                    addParams( options: [:]                                          )
include { MAFFT_ALIGN_QUERY                                      } from '../modules/local/process/mafft_align_query'                 addParams( options: modules['mafft_align_query']                 )
include { IQTREE_QUERY                                           } from '../modules/local/process/iqtree_query'                      addParams( options: modules['iqtree']                            )
include { BIOPYTHON_NEWICK_TO_NEXUS                              } from '../modules/local/process/newick_to_nexus'                   addParams( options: modules['newick_to_nexus']                   )
include { DENDROPY_LINEAGE                                       } from '../modules/local/process/dendropy_lineage'                  addParams( options: modules['assign_lineages']                   )
include { SNPS_TO_CSV                                            } from '../modules/local/process/snps_to_csv'                       addParams( options: modules['snps']                              )
include { ORDER_BY_TIPLABELS                                     } from '../modules/local/process/order_by_tiplabels'                addParams( options: modules['mafft_align_query']                 )
include { PLOT_TREE_SNPS                                         } from '../modules/local/process/plot_tree_snps'                    addParams( options: modules['plots']                             )
include { PLOT_TREE_MSA                                          } from '../modules/local/process/plot_tree_msa'                     addParams( options: modules['plots']                             )
include { REPORT as REPORT_WITH_LINEAGE                          } from '../modules/local/process/report'                            addParams( options: modules['report']                            )


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////


workflow RVFVTYPING {


    // input channels
    ch_lineages_csv              = readIn(lineages)
    ch_reps_alignment            = fasta


    // channel to collect software versions
    ch_software_versions = Channel.empty()

    
    // MODULE: Filter input files on percentage Ns
    FILTER_INPUT(ch_fasta_files)
    ch_filter_input = FILTER_INPUT.out.txt
    ch_filtered_fasta = ch_fasta_files.join(ch_filter_input)

    // remove samples that failed percentage Ns threshold
    println n_threshold
    ch_filtered_fasta
        .filter { row -> check_input( row[0], row[2], n_threshold )  }
        .map {it[0..1]}
        .set { ch_filtered_fasta }

    // MODULE: blast sequences against diamond database 

    if (params.skip_diamond) {
        ch_diamond_blastx_txt = null
    } else {
        GUNZIP_DATABASE(ch_database)
        DIAMOND_BLASTX(ch_filtered_fasta, GUNZIP_DATABASE.out.database.collect())
        ch_diamond_blastx_txt = DIAMOND_BLASTX.out.txt
        ch_diamond_blastx_txt
            .filter { row -> file(row[1])}
            .map { it[1] }
            .set { ch_diamond_blastx_txt }
        ch_software_versions = ch_software_versions.mix(DIAMOND_BLASTX.out.version.first().ifEmpty(null))
    }
    
    // MODULE: multiple sequence alignment - add query sequence to the representative sequences alignment
    MAFFT_ALIGN_QUERY (ch_filtered_fasta, ch_reps_alignment)
    ch_align_query_alignment = MAFFT_ALIGN_QUERY.out.alignment

    // MODULE: phylogenetic tree inference    
    IQTREE_QUERY (ch_align_query_alignment)
    ch_query_phylogeny = IQTREE_QUERY.out.phylogeny

    // MODULE: convert tree file from newick to nexus
    BIOPYTHON_NEWICK_TO_NEXUS (ch_query_phylogeny)
    ch_newick_to_nexus_tree = BIOPYTHON_NEWICK_TO_NEXUS.out.out_tree
    ch_software_versions = ch_software_versions.mix(BIOPYTHON_NEWICK_TO_NEXUS.out.version.first().ifEmpty(null))

    // MODULE: assign lineage
    DENDROPY_LINEAGE (ch_newick_to_nexus_tree)
    ch_software_versions = ch_software_versions.mix(DENDROPY_LINEAGE.out.version.first().ifEmpty(null))

    // MODULE: convert masked snp alignment to csv
    SNPS_TO_CSV ( ch_align_query_alignment )
    ch_snps_to_csv = SNPS_TO_CSV.out.csv

    // MODULE: order alignment by tip labels and plot tree with msa
    ch_tree_alignment = ch_newick_to_nexus_tree.join(ch_align_query_alignment)
    ch_tree_alignment
        .map { row -> 
        def meta = [:]
        meta.id = row[0].id
        [ meta , [ file(row[1], checkIfExists: true), file(row[2], checkIfExists: true) ] ]
    }.set { ch_tree_alignment }

    ORDER_BY_TIPLABELS ( ch_tree_alignment )

    ch_order_by_tiplabels_align_fasta = ORDER_BY_TIPLABELS.out.fasta

    // get lineages channel for plotting
    ch_metadata = ch_lineages_csv.map{row -> row[1]}
    
    // MODULE: plot trees with snps data, lineages metadata and multiple sequence alignment
    ch_tree_snps_msa_input = ch_query_phylogeny.join(ch_snps_to_csv).join(ch_order_by_tiplabels_align_fasta)
    ch_tree_snps_msa_input
        .map { row -> 
        def meta = [:]
        meta.id = row[0].id
        [ meta , [ file(row[1], checkIfExists: true), file(row[2], checkIfExists: true), file(row[3], checkIfExists: true) ] ]
    }.set { ch_tree_snps_msa_input }


    // MODULE: Plot
    PLOT_TREE_SNPS ( ch_tree_snps_msa_input, ch_metadata.collect())
    PLOT_TREE_MSA ( ch_tree_snps_msa_input )


    // MODULE: Report results
    ch_filtered_fasta
        .filter { row -> row[1]}
        .map { it[1] }
        .set { ch_filtered_fasta }

    ch_assign_lineage_csv = DENDROPY_LINEAGE.out.lineage_csv 
    ch_assign_lineage_csv
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_assign_lineage_csv }
    
    REPORT_WITH_LINEAGE (
        ch_filtered_fasta.collect(), 
        ch_assign_lineage_csv.collect(), 
        ch_lineages_csv, 
        ch_diamond_blastx_txt.collect()
        )


    /*
     * MODULE: Pipeline reporting
     */
    
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    // Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/