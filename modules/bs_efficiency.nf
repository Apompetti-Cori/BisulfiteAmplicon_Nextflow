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
params.pubdir = "bs_efficiency"

process bs_efficiency {
    maxForks 3
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id),
    path(chg_ot),
    path(chg_ob),
    path(chh_ot),
    path(chh_ob)

    output:
    tuple val(file_id), path("*${file_id}*.tsv")

    shell:
    '''
    #!/usr/bin/env Rscript

    library(magrittr)

    ### read in files containing location of non-CpG methylation
    #CHG context
    rbind(vroom::vroom("!{chg_ot}", delim = '\t', skip = 1, 
                   col_names = c('read_id', 'strand', 'chr', 'position', 
                                 'meth_status')),
            vroom::vroom("!{chg_ob}", delim = '\t', skip = 1, 
                   col_names = c('read_id', 'strand', 'chr', 'position', 
                                 'meth_status'))) %>%
    dplyr::mutate(context = 'CHG') -> chg
    #CHH context
    rbind(vroom::vroom("!{chh_ot}", delim = '\t', skip = 1, 
                   col_names = c('read_id', 'strand', 'chr', 'position', 
                                 'meth_status')),
            vroom::vroom("!{chh_ob}", delim = '\t', skip = 1, 
                   col_names = c('read_id', 'strand', 'chr', 'position', 
                                 'meth_status'))) %>%
    dplyr::mutate(context = 'CHH') -> chh

    ### count percentage of "methylated" non-CpG cytosines
    rbind(chg, chh) %>%
    #distinct() %>%
    dplyr::mutate(meth_status = ifelse(meth_status %in% c('x', 'h'), 
                                     'unmeth', 'meth')) %>%
    dplyr::count(context, meth_status) %>%
    tidyr::pivot_wider(names_from = meth_status, values_from = n) %>%
    tidyr::replace_na(list(meth = 0, unmeth = 0)) %>%
    dplyr::mutate(percent_meth = (meth / (meth + unmeth)) * 100) -> bs_efficiency

    #save
    readr::write_tsv(bs_efficiency, "!{file_id}_bs_efficiency.tsv")
    '''
}