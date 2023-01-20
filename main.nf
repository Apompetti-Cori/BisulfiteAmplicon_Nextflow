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
params.genome = false
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false

//Include modules to main pipeline
include { fastqc as pretrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'pretrim_fastqc')
include { trim_galore } from './modules/trim_galore.nf'
include { fastqc as posttrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'posttrim_fastqc')
include { bismark_align } from './modules/bismark_align.nf' addParams(db: params.db)
include { bismark_extract } from './modules/bismark_extract.nf'

//Create channel for reads. By default, auto-detects paired end data. Specify --singleEnd if your fastq files are in single-end format
Channel
.fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2)
.ifEmpty {exit 1, "Cannot find any reads matching ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line."}
.set{ reads_ch }

workflow {

    //Run fastqc on raw reads
    pretrim_fastqc(reads_ch)
    //Run trim_galore on raw reads
    trim_galore(reads_ch)

    //Run fastqc on trimmed reads, specifies trim_galore[0] because second input channel is not need for this process
    posttrim_fastqc(trim_galore.out[0])
    //Run bismark_align on trimmed reads
    bismark_align(trim_galore.out[0])

    //Run bismark_extract on bismark_align output
    bismark_extract(bismark_align.out)
}

