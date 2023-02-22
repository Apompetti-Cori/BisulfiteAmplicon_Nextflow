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
params.pubdir = "conv_stats/merge"
params.lambda_rname = "J02459.1_lambda"

process CONV_STATS_MERGE {
    memory '8 GB'
    cpus 1

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    path('conv_stats/*')

    output:
    path('conv_stats_merge.tsv')

    shell:
    '''
    #!/bin/bash

    awk '(NR == 1) || (FNR > 1)' *.tsv > conv_stats_merge.tsv
    '''
}