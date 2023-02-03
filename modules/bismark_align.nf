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
params.pubdir = "bismark_align"
params.db = false

process bismark_align {
    maxForks 10
    memory '32 GB'
    cpus 4
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(reads)

    output:
    tuple val(file_id), path("*.bam")
    path("*report.txt"), emit: align_report

    script:
    """
    bismark --bowtie2 --parallel ${task.cpus} ${params.db} -1 ${reads[0]} -2 ${reads[1]}
    """
}