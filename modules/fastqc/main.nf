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
process FASTQC {
    maxForks 4
    memory '8 GB'
    cpus 2
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.batch == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'link'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.{html,zip}")

    script:
    """
    fastqc -t $task.cpus $reads
    """
}