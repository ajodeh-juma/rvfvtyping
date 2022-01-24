////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// def valid_params = [
//     group            : ['viral', 'bacteria', 'archaea', 'fungi'],
//     formats          : ['fasta', 'protein-fasta'],
//     section          : ['refseq', 'genbank'],
//     assembly_levels  : ['all', 'complete', 'chromosome', 'scaffold', 'contig'],
//     guide_tree       : ['ml-snps', 'ml-bayes']
// ]

params.summary_params = [:]

// Validate input parameters
// def groups = params.group ? params.group.split(',').collect{ it.trim().toLowerCase() } : []
// if (!valid_params['group'].contains(params.group)) {
//     exit 1, "Invalid group option: ${params.group}. Valid options for '--group': ${valid_params['group'].join(', ')}."
// }

// def formats = params.formats ? params.formats.split(',').collect{ it.trim().toLowerCase() } : []
// if (!valid_params['formats'].contains(params.formats)) {
//     exit 1, "Invalid formats option: ${params.formats}. Valid options for '--formats': ${valid_params['formats'].join(', ')}."
// }

// def section = params.section ? params.section.split(',').collect{ it.trim().toLowerCase() } : []
// if (!valid_params['section'].contains(params.section)) {
//     exit 1, "Invalid section option: ${params.section}. Valid options for '--section': ${valid_params['section'].join(', ')}."
// }

// def assembly_levels = params.assembly_levels ? params.assembly_levels.split(',').collect{ it.trim().toLowerCase() } : []
// if (!valid_params['assembly_levels'].contains(params.assembly_levels)) {
//     exit 1, "Invalid assembly_levels option: ${params.section}. Valid options for '--assembly_levels': ${valid_params['assembly_levels'].join(', ')}."
// }

// def guide_treesList = ['ml-snps', 'ml-bayes']
// if (!guide_treesList.contains(params.guide_tree)) {
//     exit 1, "Invalid guide tree algorithm option: ${params.guide_tree}. Valid options: ${guide_treesList.join(', ')}"
// }


// Check input path parameters to see if they exist
checkPathParamList = [
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// // Check mandatory parameters
// if (!params.download_genomes && !params.db) { exit 1, "No download RefSeq option or database specified!" }
// if (params.download_genomes && params.db)   { Checks.download_genomes_db_warn(log) }
// if (params.download_genomes) {Checks.download_genomes_warn(log)}
// if (params.db) { ch_database = file(params.db) } 

// // group channel
// if (params.group) {
//     ch_group = Channel
//         .from(params.group)
// }

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

// urls channels
// ch_taxonomy_url = Channel.from(params.taxonomy_url)
// ch_protein_acc2taxid_url = Channel.from(params.protein_accession2taxid_url)



// Stage required files

// known lineages file
lineages_csv  = file("$projectDir/assets/G2-M-segment_Lineages.csv", checkIfExists: true)

// diamond database file
// db = file("$projectDir/db/viral.protein.faa.dmnd.gz")
ch_database = Channel.fromPath(params.db)
    .map { fn -> 
           def meta = [:]
           meta.id = fn.baseName
           path = fn.parent
           db = [ meta, fn] 
           return  db 
         }

// representative sequences alignments
reps_aln = file("$projectDir/assets/rvfv.representative.fasta", checkIfExists: true)

// representative sequences guide trees
// reps_tree = file("$projectDir/assets/rvfv.representative.fasta.treefile", checkIfExists: true)

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
// include { SUMMARY_LINEAGES as REPORT_WITH_LINEAGE                } from '../modules/local/process/summarize_lineages'                addParams( options: modules['summary']                           )

include { REPORT as REPORT_WITH_LINEAGE                          } from '../modules/local/process/report'                            addParams( options: modules['report']                            )


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////


workflow RVFVTYPING {


    // lineages file channel
    ch_lineages_csv                  = readIn(lineages_csv)

    // channels for the different reference representative files
    ch_reps_alignment                = reps_aln

    // channel to collect software versions
    ch_software_versions = Channel.empty()

    
    // MODULE: Filter input files on percentage Ns
    FILTER_INPUT(ch_fasta_files)
    ch_filter_input = FILTER_INPUT.out.txt
    ch_filtered_fasta = ch_fasta_files.join(ch_filter_input)

    // remove samples that failed percentage Ns threshold
    ch_filtered_fasta
        .filter { row -> check_input( row[0], row[2], 50.0 )  }
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
    ch_order_by_tiplabels_input = ch_newick_to_nexus_tree.join(ch_align_query_alignment)
    ch_order_seq_tree = ch_order_by_tiplabels_input
        .map { row -> 
        def meta = [:]
        tree = row[1]
        meta.id = row[0].id
        [ meta, tree ]
        }
    
    ch_order_seq_align = ch_order_by_tiplabels_input
        .map { row -> 
        def meta = [:]
        align = row[2]
        meta.id = row[0].id
        [ meta, align ]
        }

    ORDER_BY_TIPLABELS ( ch_order_seq_tree, ch_order_seq_align )
    ch_order_by_tiplabels_align_fasta = ORDER_BY_TIPLABELS.out.fasta

    // get lineages channel for plotting
    ch_metadata = ch_lineages_csv.map{row -> row[1]}
    
    // MODULE: plot trees with snps data, lineages metadata and multiple sequence alignment
    ch_tree_snps_msa_input = ch_query_phylogeny.join(ch_snps_to_csv).join(ch_order_by_tiplabels_align_fasta)
    ch_tree = ch_tree_snps_msa_input
        .map { row -> 
        def meta = [:]
        tree = row[1]
        meta.id = row[0].id
        [ meta, tree ]
        }
    ch_snps = ch_tree_snps_msa_input
        .map { row -> 
        def meta = [:]
        snps = row[2]
        meta.id = row[0].id
        [ meta, snps  ]
        }
    ch_msa = ch_tree_snps_msa_input
        .map { row -> 
        def meta = [:]
        msa = row[3]
        meta.id = row[0].id
        [ meta, msa  ]
        }
        
    PLOT_TREE_SNPS ( ch_tree, ch_snps, ch_metadata.collect())
    PLOT_TREE_MSA ( ch_tree, ch_msa ) 


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