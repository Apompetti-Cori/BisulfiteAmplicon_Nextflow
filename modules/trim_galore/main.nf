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

/*
Run trim_galore on each read stored within the reads_ch channel
*/
process TRIM_GALORE {
    maxForks 4
    memory '8 GB'
    cpus 4

    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.batch == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'link'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path("*trimming_report.txt"), emit: report

    script:
    def singleEnd = meta.single_end ? '' : '--paired'
    def rrbs = params.rrbs ? '--rrbs' : ''

    """
    trim_galore \
    ${singleEnd} \
    ${rrbs} \
    --length 35 \
    --quality 28 \
    --phred33 \
    --cores ${task.cpus} \
    ${reads}
    """
}