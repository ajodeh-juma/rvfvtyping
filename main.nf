#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rvfvtyping
========================================================================================
 rvfvtyping Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/ajodeh-juma/rvfvtyping
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run main.nf --input '*.fasta' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)


////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    include { RVFVTYPING } from './workflows/rvfvtyping' addParams( summary_params: summary_params )
    RVFVTYPING ()
}
