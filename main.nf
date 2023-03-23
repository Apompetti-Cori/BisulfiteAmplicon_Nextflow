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
params.fastq_folder = "${workflow.projectDir}/fastq"
params.reads = "${params.fastq_folder}/*{_L00}{1,2,3,4}{_R,_}{1,2}*.{fastq,fq}.gz"
params.singleEnd = false
params.multiLane = false
params.rrbs = false
params.amplifytargets = false
params.genome = false
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false
params.cpg_wl = "${workflow.projectDir}/tables/3_cpg_whitelist.tsv"
params.ref_dist = "${workflow.projectDir}/tables/2021-07_new_target_cb_ref.tsv"

//Include modules to main pipeline
include { FASTQC as PRETRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'pretrim_fastqc')
include { TRIM_GALORE } from './modules/trim_galore/main'
include { FASTQC as POSTTRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'posttrim_fastqc')
include { BISMARK_ALIGN } from './modules/bismark/bismark_align/main' addParams(db: params.db)
include { BISMARK_EXTRACT } from './modules/bismark/bismark_extract/main'
include { CONV_STATS_CREATE } from './modules/conversion_stats/create_stats/main'
include { MULTIQC } from './modules/multiqc/main'
include { bs_efficiency } from './modules/bs_efficiency.nf'
include { allele_freq } from './modules/allele_freq.nf'
include { calc_summary } from './modules/calc_summary.nf'

//Provide sample table in csv format to have pipeline process samples via sample table
params.sample_table = false
if ( params.sample_table ){
    // Channel for the samplesheet
    ch_samplesheet = Channel.fromPath(params.sample_table)

    // Parse it line by line
    reads_ch = ch_samplesheet.splitCsv(header:true).map {

        // This is the read1 and read2 entry
        r1_L1 = it['r1_L1']
        r1_L2 = it['r1_L2']
        r1_L3 = it['r1_L3']
        r1_L4 = it['r1_L4']
        r2_L1 = it['r2_L1']
        r2_L2 = it['r2_L2']
        r2_L3 = it['r2_L3']
        r2_L4 = it['r2_L4']

        // Detect wiether single-end or paired-end
        is_singleEnd = r2_L1.toString()=='' ? true : false
        is_multiLane = r1_L2.toString()=='' ? false : true
        
        // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
        meta = [id: it['sample'], single_end: is_singleEnd, multi_lane: is_multiLane]
        
        // We return a nested map, the first entry is the meta map, the second one is the read(s)
        if ( is_singleEnd ){
            if ( r1_L4.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], []]
            }
            else if ( r1_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], []]
            }
            else if ( r1_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], []]
            }
            else {
                [meta, [r1_L1], []]
            }
        }
        else{
            if (r2_L4.toString()!=''){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], [r2_L1, r2_L2, r2_L3, r2_L4]]
            }
            else if ( r2_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], [r2_L1, r2_L2, r2_L3]]
            }
            else if ( r2_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], [r2_L1, r2_L2]]
            }
            else {
                [meta, [r1_L1], [r2_L1]]
            }
        }
    }
}
else {
    Channel
    .fromFilePairs(params.reads, size: params.singleEnd ? (params.multiLane ?: 1) : 2)
    .ifEmpty {exit 1, "Cannot find any reads matching ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line."}
    .set{ reads_ch }
}


process LANE_COMBINE {
    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*{_mL,_sL,_sL_1,_sL_2,_mL_1,_mL_2}.fq.gz")

    script:
    if ( meta.multi_lane ){
        if (meta.single_end){
            """
            cat ${reads_1} > ${meta.id}_mL.fq.gz
            """
        }
        else {
            """
            cat ${reads_1} > ${meta.id}_mL_1.fq.gz
            cat ${reads_2} > ${meta.id}_mL_2.fq.gz
            """
        }
    }
    else {
        if (meta.single_end){
            """
            cat ${reads_1} > ${meta.id}_sL.fq.gz
            """
        }
        else {
            """
            cat ${reads_1} > ${meta.id}_sL_1.fq.gz
            cat ${reads_2} > ${meta.id}_sL_2.fq.gz
            """
        }
    }
}


workflow {
    LANE_COMBINE(reads_ch)
    reads_ch = LANE_COMBINE.out
    reads_ch.view()

    
    //Run fastqc on raw reads
    PRETRIM_FASTQC(reads_ch)
    //Run trim_galore on raw reads
    TRIM_GALORE(reads_ch)

    //Run fastqc on trimmed reads, specifies trim_galore[0] because second input channel is not need for this process
    POSTTRIM_FASTQC(TRIM_GALORE.out.reads)
    //Run bismark_align on trimmed reads
    //State Dependency: Wait until TRIM_GALORE is done to run bismark align
    //State Dependency: Wait until POSTTRIM_FASTQC is done to run bismark align
    state = POSTTRIM_FASTQC.out.collect()
    BISMARK_ALIGN(TRIM_GALORE.out.reads.collect(flat: false).flatMap(), state)

    //Run bismark_extract on bismark_align output
    BISMARK_EXTRACT(BISMARK_ALIGN.out.bam.collect(flat: false).flatMap())

    //Run bisulfite_conversion on bismark_align output
    CONV_STATS_CREATE(BISMARK_ALIGN.out.bam.collect(flat: false).flatMap())

    //Run multiqc on pretrim fastqc output, trim_galore trimming report, posttrim fastqc output, bismark conversion output
    MULTIQC(PRETRIM_FASTQC.out.collect().combine(POSTTRIM_FASTQC.out.collect()).combine(TRIM_GALORE.out.report.collect()).combine(BISMARK_ALIGN.out.report.collect()).combine(CONV_STATS_CREATE.out.report.collect()))
    
    //Run bs_efficiency on bismark_extract chg (ot,ob) and chh (ot,ob) output
    //bs_efficiency(bismark_extract.out.chg_ot.combine(bismark_extract.out.chg_ob, by: 0).combine(bismark_extract.out.chh_ot.combine(bismark_extract.out.chh_ob, by: 0), by: 0))

    if( params.amplifytargets ){
        //Run allele_freq on bismark_extract cpg (ot,ob) output
        allele_freq(bismark_extract.out.cpg_ot.combine(bismark_extract.out.cpg_ob, by: 0).combine(Channel.fromPath( "${params.cpg_wl}" )))

        //Run calc_summary on allele_freq and bs_efficiency output
        calc_summary(bs_efficiency.out.combine(allele_freq.out, by: 0).combine(Channel.fromPath( "${params.ref_dist}" )))
    }
}

