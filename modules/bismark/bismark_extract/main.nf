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

process BISMARK_EXTRACT {
    maxForks 4
    memory '16 GB'
    cpus 4

    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.batch == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'link'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bismark.cov.gz"), optional: true
    tuple val(meta), path("CpG_OT_${meta.id}*"), emit: cpg_ot, optional: true
    tuple val(meta), path("CpG_OB_${meta.id}*"), emit: cpg_ob, optional: true
    tuple val(meta), path("CHG_OT_${meta.id}*"), emit: chg_ot, optional: true
    tuple val(meta), path("CHG_OB_${meta.id}*"), emit: chg_ob, optional: true
    tuple val(meta), path("CHH_OT_${meta.id}*"), emit: chh_ot, optional: true
    tuple val(meta), path("CHH_OB_${meta.id}*"), emit: chh_ob, optional: true

    script:
    if( meta.single_end ){
        """
        bismark_methylation_extractor --parallel ${task.cpus} --gzip --bedGraph ${bam}
        """
    }
    else{
        """
        bismark_methylation_extractor --parallel ${task.cpus} --paired-end --include_overlap --gzip --bedGraph ${bam}
        """
    }

}
