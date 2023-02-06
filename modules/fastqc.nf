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
params.pubdir = "fastqc"

/*
Run fastqc on fastq files
*/
process fastqc {
    maxForks 10
    memory '8 GB'
    cpus 2
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(reads)

    output:
    path("*.{html,zip}")

    script:
    """
    fastqc -t $task.cpus $reads
    """
}