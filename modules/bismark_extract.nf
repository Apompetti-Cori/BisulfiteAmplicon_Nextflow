#!/usr/bin/env nextflow

/*
Coriell Institute for Medical Research
Bisulfite Amplicon Pipeline. Started January 2023.

Contributors:
Anthony Pompetti <apompetti@coriell.org>

Methodology adapted from:
prior snakemake pipeline developed by Matthew Walt
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "bismark_extract"

process bismark_extract {
    maxForks 3
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(reads)

    output:
    tuple val(file_id), path("*${file_id}*")

    script:
    """
    bismark_methylation_extractor --multicore ${task.cpus} --paired-end --include_overlap --bedGraph ${reads}
    """
}
