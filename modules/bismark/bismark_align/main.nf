#!/usr/bin/env nextflow

/*
================================================================================
Coriell Institute for Medical Research

Contributors:
Anthony Pompetti <apompetti@coriell.org>
================================================================================
*/

/*
================================================================================
Enable Nextflow DSL2
================================================================================
*/
nextflow.enable.dsl=2

/*
================================================================================
Configurable variables for module
================================================================================
*/
params.outdir = "./results"
params.pubdir = "bismark_align"
params.db = false

/*
================================================================================
Module declaration
================================================================================
*/
process BISMARK_ALIGN {
    maxForks 4
    memory '40 GB'
    cpus 4
    
    // Check batch and save output accordingly
    publishDir "${params.outdir}",  saveAs: { meta.batch == '' ? "${params.pubdir}/${it}" : "${meta.batch}/${params.pubdir}/${it}" }, mode: 'link'

    input:
    tuple val(meta), path(reads)
    val(state)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path("*report.txt"), emit: report

    script:
    if(meta.single_end){
        """
        bismark --bowtie2 --parallel ${task.cpus} ${params.db} ${reads}
        """
    }
    else{
        """
        bismark --bowtie2 --parallel ${task.cpus} ${params.db} -1 ${reads[0]} -2 ${reads[1]}
        """
    }

}