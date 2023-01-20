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

//Configurable variables for pipeline
params.reads = "${workflow.projectDir}/fastq/*{_R,_}{1,2}*.{fastq,fq}.gz"
params.singleEnd = false
params.multiqc_config = "${workflow.projectDir}/multiqc_config.yaml"

//Include modules to main pipeline
include { fastqc as pretrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'pretrim_fastqc')
include { trim_galore } from './modules/trim_galore.nf'

//Create channel for reads. By default, auto-detects paired end data. Specify --singleEnd if your fastq files are in single-end format
Channel
.fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2)
.ifEmpty {exit 1, "Cannot find any reads matching ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line."}
.set{ reads_ch }

workflow {

    //Perform fastqc on raw reads, trim the reads with trim_galore, and perform fastqc on trimmed reads
    pretrim_fastqc(reads_ch)
    trim_galore(reads_ch)

}

