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
params.pubdir = "allele_freq"

process allele_freq {
    maxForks 3
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id),
    path(cpg_ot),
    path(cpg_ob),
    path(cpg_wl)

    output:
    tuple val(file_id), path("${file_id}*.tsv")

    shell:
    '''
    #!/usr/bin/env Rscript
    library(magrittr)
    library(tidyverse)
    library(vroom)

    ### Read in files and whitelisted CpGs
    vroom("!{cpg_ot}", delim = '\t', skip = 1, 
        col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')) %>%
    select(-strand) %>%
    distinct() -> cpg_ot
    
    vroom("!{cpg_ob}", 
        delim = '\t', skip = 1,
        col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')) %>%
    select(-strand) %>%
    distinct() -> cpg_ob

    vroom("!{cpg_wl}") -> cpg_whitelist

    ### Count the number of methylated CpGs per read
    rbind(cpg_ot, cpg_ob) %>%
        unite(chr_base, c(chr, position), sep = '_') %>%
        left_join(cpg_whitelist, by = 'chr_base') %>%
        na.omit() %>%
        add_count(read_id, target) %>%
        filter(n == target_cpg_count) %>%
        select(-n) %>%
        count(read_id, target, target_cpg_count, meth_status) %>%
        pivot_wider(names_from = meth_status, values_from = n) %>%
        replace_na(list(z = 0, Z = 0)) %>%
        select(read_id, target, target_cpg_count, 
            count_meth_cpgs = Z) -> read_counts
    #write_tsv(read_counts, "!{file_id}_summarized_counts.tsv")

    ### Summarize counts by target and number of methylated CpGs
    tbl <- NULL
    for (i in unique(cpg_whitelist$target)) {
    rbind(tibble(target = i,
        count_meth_cpgs = 0:unique(filter(cpg_whitelist, 
        target == i)$target_cpg_count)), 
        tbl) -> tbl
    }

    read_counts %>%
        add_count(target, name = 'target_reads') %>% 
        count(target, target_reads, count_meth_cpgs, 
            name = 'read_count') %>%
        mutate(probability = read_count / target_reads) %>%
        full_join(tbl, by = c('target', 'count_meth_cpgs')) %>%
        arrange(target, count_meth_cpgs) %>%
        fill(target_reads) %>%
        replace_na(list(read_count = 0, probability = 0)) -> counts
    write_tsv(counts, "!{file_id}_summarized_counts.tsv")
    '''
}