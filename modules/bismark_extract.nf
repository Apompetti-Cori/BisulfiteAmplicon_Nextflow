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
    maxForks 10
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(bam)

    output:
    tuple val(file_id), path("*.bismark.cov.gz"), optional: true
    tuple val(file_id), path("CpG_OT_${file_id}*"), emit: cpg_ot, optional: true
    tuple val(file_id), path("CpG_OB_${file_id}*"), emit: cpg_ob, optional: true
    tuple val(file_id), path("CHG_OT_${file_id}*"), emit: chg_ot, optional: true
    tuple val(file_id), path("CHG_OB_${file_id}*"), emit: chg_ob, optional: true
    tuple val(file_id), path("CHH_OT_${file_id}*"), emit: chh_ot, optional: true
    tuple val(file_id), path("CHH_OB_${file_id}*"), emit: chh_ob, optional: true

    script:
    """
    bismark_methylation_extractor --multicore ${task.cpus} --paired-end --include_overlap --bedGraph ${bam}
    """
}
