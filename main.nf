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
params.reads = "${params.fastq_folder}/*{_R,_}{1,2}*.{fastq,fq}.gz"
params.singleEnd = false
params.rrbs = false
params.amplifytargets = false
params.genome = false
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false
params.cpg_wl = "${workflow.projectDir}/tables/3_cpg_whitelist.tsv"
params.ref_dist = "${workflow.projectDir}/tables/2021-07_new_target_cb_ref.tsv"

//Include modules to main pipeline
include { fastqc as pretrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'pretrim_fastqc')
include { trim_galore } from './modules/trim_galore.nf'
include { fastqc as posttrim_fastqc } from './modules/fastqc.nf' addParams(pubdir: 'posttrim_fastqc')
include { bismark_align } from './modules/bismark_align.nf' addParams(db: params.db)
include { bismark_extract } from './modules/bismark_extract.nf'
include { bisulfite_conversion } from './modules/bisulfite_conversion.nf'
include { bs_efficiency } from './modules/bs_efficiency.nf'
include { allele_freq } from './modules/allele_freq.nf'
include { calc_summary } from './modules/calc_summary.nf'
include { multiqc } from './modules/multiqc.nf'

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

    //Run bisulfite_conversion on bismark_align output
    bisulfite_conversion(bismark_align.out)

    //Run multiqc on pretrim fastqc output, trim_galore trimming report, posttrim fastqc output, bismark conversion output
    multiqc(pretrim_fastqc.out.collect().combine(posttrim_fastqc.out.collect()).combine(trim_galore.out.trimming_report.collect()).combine(bisulfite_conversion.out.collect()))

    //Run bs_efficiency on bismark_extract chg (ot,ob) and chh (ot,ob) output
    bs_efficiency(bismark_extract.out.chg_ot.combine(bismark_extract.out.chg_ob, by: 0).combine(bismark_extract.out.chh_ot.combine(bismark_extract.out.chh_ob, by: 0), by: 0))

    if( params.amplifytargets ){
        //Run allele_freq on bismark_extract cpg (ot,ob) output
        allele_freq(bismark_extract.out.cpg_ot.combine(bismark_extract.out.cpg_ob, by: 0).combine(Channel.fromPath( "${params.cpg_wl}" )))

        //Run calc_summary on allele_freq and bs_efficiency output
        calc_summary(bs_efficiency.out.combine(allele_freq.out, by: 0).combine(Channel.fromPath( "${params.ref_dist}" )))
    }
}

