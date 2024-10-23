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
Configurable variables for pipeline
================================================================================
*/
params.rrbs = false
params.input_type = "fastq"
params.genome = false
params.sample_table = false //Provide sample table in csv format to have pipeline process samples via sample table
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false

/*
================================================================================
Include functions to main pipeline
================================================================================
*/
include { createInputChannel } from './modules/preprocess/functions.nf'

/*
================================================================================
Include modules to main pipeline
================================================================================
*/
include { PREPROCESS_READS as PREPROCESS } from './modules/preprocess/main.nf'
include { FASTQC as PRETRIM_FASTQC } from './modules/fastqc/main.nf' addParams(pubdir: 'pretrim_fastqc')
include { TRIM_GALORE } from './modules/trim_galore/main.nf'
include { FASTQC as POSTTRIM_FASTQC } from './modules/fastqc/main.nf' addParams(pubdir: 'posttrim_fastqc')
include { BISMARK_ALIGN } from './modules/bismark/bismark_align/main.nf' addParams(db: params.db)
include { BISMARK_EXTRACT } from './modules/bismark/bismark_extract/main.nf'
include { CONV_STATS_CREATE } from './modules/conversion_stats/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

/*
================================================================================
Channel creation
================================================================================
*/
index_ch = Channel.fromPath(
  params.db, type: 'dir'
)

/*
================================================================================
Workflow declaration
================================================================================
*/
workflow {
    input_ch = createInputChannel(params.sample_table, params.input_type)
    PREPROCESS(reads_ch)
    reads_ch = PREPROCESS.out

    index_ch.view()

    //Run fastqc on raw reads
    PRETRIM_FASTQC(reads_ch)
    
    //Run trim_galore on raw reads
    TRIM_GALORE(reads_ch)

    //Run fastqc on trimmed reads, specifies trim_galore[0] because second input channel is not need for this process
    POSTTRIM_FASTQC(TRIM_GALORE.out.reads)

    //State Dependency: Wait until POSTTRIM_FASTQC is done to run bismark align
    state = POSTTRIM_FASTQC.out.collect()

    //Run bismark_align on trimmed reads
    BISMARK_ALIGN(TRIM_GALORE.out.reads.collect(flat: false).flatMap(), state)

    //Run bismark_extract on bismark_align output
    BISMARK_EXTRACT(BISMARK_ALIGN.out.bam.collect(flat: false).flatMap())

    //Run bisulfite_conversion on bismark_align output
    CONV_STATS_CREATE(BISMARK_ALIGN.out.bam.collect(flat: false).flatMap())

    //Run multiqc on pretrim fastqc output, trim_galore trimming report, posttrim fastqc output, bismark conversion output
    MULTIQC(PRETRIM_FASTQC.out.collect().combine(POSTTRIM_FASTQC.out.collect()).combine(TRIM_GALORE.out.report.collect()).combine(BISMARK_ALIGN.out.report.collect()).combine(CONV_STATS_CREATE.out.report.collect()))
}

