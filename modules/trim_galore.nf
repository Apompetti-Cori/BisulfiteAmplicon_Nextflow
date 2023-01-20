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
params.pubdir = "trim_galore"
params.singleEnd = false

/*
Run trim_galore on each read stored within the reads_ch channel
*/
process trim_galore {
    maxForks 3
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    tuple val(file_id), path(reads)

    output:
    tuple val(file_id), path("*.fq.gz")
    path("*trimming_report.txt"), emit: trimming_report

    script:
    if ( params.singleEnd )
    """
    trim_galore \
    --length 35 \
    --quality 28 \
    --phred33 \
    --cores ${task.cpus} \
    $reads
    """

    else
    """
    trim_galore \
    --length 35 \
    --quality 28 \
    --paired \
    --phred33 \
    --clip_R2 3 \
    --cores ${task.cpus} \
    $reads
    """
}