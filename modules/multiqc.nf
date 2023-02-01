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
params.pubdir = "multiqc"
params.multiqc_fn = "multiqc_report"

process multiqc {
    memory '32 GB'
    cpus 8

    conda '/opt/miniconda3/envs/multiqc'

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'
    
    input:
    path(multiqc_config)
    path('multiqc_input/*')

    output:
    path("${params.multiqc_fn}.html")

    script:
    """
    multiqc multiqc_input/ --config ${multiqc_config} --filename ${params.multiqc_fn}
    """
}