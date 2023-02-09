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
params.pubdir = "bs_conversion"
params.lambda_rname = "J02459.1_lambda"

process bisulfite_conversion {
    maxForks 4
    memory '8 GB'
    cpus 1

    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(file_id), path(bam)

    output:
    tuple val(file_id), path("*.tsv"), optional: true
    path("*.tsv"), emit: conversion_report, optional: true

    shell:
    $/
    #!/usr/bin/env python

    from pathlib import Path
    from collections import Counter
    import pysam
    import re

    if Path("!{bam}").suffix == ".bam":
        bamfile = pysam.AlignmentFile(str("!{bam}"), "rb")
    elif Path("!{bam}").suffix == ".sam":
        bamfile = pysam.AlignmentFile(str("!{bam}"), "r")
    else:
        raise ValueError("Alignment file must end in .sam or .bam")

    sample_id = re.sub('_val.*', '', str("!{file_id}"))
    
    outfile = str(sample_id) + str(".conversion-stats.tsv")

    counts = {"z":0, "Z":0, "x":0, "X":0, "H":0, "h":0, "H":0, "u":0, "U":0}
    lambda_counts =  {"z":0, "Z":0, "x":0, "X":0, "H":0, "h":0, "H":0, "u":0, "U":0}

    for read in bamfile:
        ref = read.reference_name
        aln_str = read.get_tag("XM")
        if not read.is_unmapped:
            if (read.is_paired and read.is_proper_pair):
                c = Counter(aln_str)
                if ref == "!{params.lambda_rname}":
                    lambda_counts["z"] += c.get("z", 0)
                    lambda_counts["Z"] += c.get("Z", 0)
                    lambda_counts["x"] += c.get("x", 0)
                    lambda_counts["X"] += c.get("X", 0)
                    lambda_counts["h"] += c.get("h", 0)
                    lambda_counts["H"] += c.get("H", 0)
                    lambda_counts["u"] += c.get("u", 0)
                    lambda_counts["U"] += c.get("U", 0)
                counts["z"] += c.get("z", 0)
                counts["Z"] += c.get("Z", 0)
                counts["x"] += c.get("x", 0)
                counts["X"] += c.get("X", 0)
                counts["h"] += c.get("h", 0)
                counts["H"] += c.get("H", 0)
                counts["u"] += c.get("u", 0)
                counts["U"] += c.get("U", 0)
            else: #Mapped single end reads
                c = Counter(aln_str)
                if ref == lambda_rname:
                    lambda_counts["z"] += c.get("z", 0)
                    lambda_counts["Z"] += c.get("Z", 0)
                    lambda_counts["x"] += c.get("x", 0)
                    lambda_counts["X"] += c.get("X", 0)
                    lambda_counts["h"] += c.get("h", 0)
                    lambda_counts["H"] += c.get("H", 0)
                    lambda_counts["u"] += c.get("u", 0)
                    lambda_counts["U"] += c.get("U", 0)
                counts["z"] += c.get("z", 0)
                counts["Z"] += c.get("Z", 0)
                counts["x"] += c.get("x", 0)
                counts["X"] += c.get("X", 0)
                counts["h"] += c.get("h", 0)
                counts["H"] += c.get("H", 0)
                counts["u"] += c.get("u", 0)
                counts["U"] += c.get("U", 0)
    bamfile.close()

    # Tabulate counts and calculate percentages
    unmeth = counts["z"] + counts["x"] + counts["h"] + counts["u"]
    meth = counts["Z"] + counts["X"] + counts["H"] + counts["U"]
    perc_meth = round(meth / (meth + unmeth) * 100, 3)
    cg_meth = counts["Z"]
    cg_unmeth = counts["z"]
    perc_cg_meth = round(cg_meth / (cg_meth + cg_unmeth) * 100, 3)
    nonCg_meth = counts["X"] + counts["H"] + counts["U"]
    nonCg_unmeth = counts["x"] + counts["h"] + counts["u"]
    perc_nonCg_meth = round(nonCg_meth / (nonCg_meth + nonCg_unmeth) * 100, 3)
    lambda_unmeth = lambda_counts["z"] + lambda_counts["x"] + lambda_counts["h"] + lambda_counts["u"]
    lambda_meth = lambda_counts["Z"] + lambda_counts["X"] + lambda_counts["H"] + lambda_counts["U"]
    perc_lambda_meth = round(lambda_meth / (lambda_meth + lambda_unmeth) * 100, 3)
    lambda_cg_meth = lambda_counts["Z"]
    lambda_cg_unmeth = lambda_counts["z"]
    perc_lambda_cg_meth = round(lambda_cg_meth / (lambda_cg_meth + lambda_cg_unmeth) * 100, 3)
    lambda_nonCg_meth = lambda_counts["X"] + lambda_counts["H"] + lambda_counts["U"]
    lambda_nonCg_unmeth = lambda_counts["x"] + lambda_counts["h"] + lambda_counts["u"]
    perc_lambda_nonCg_meth = round(lambda_nonCg_meth / (lambda_nonCg_meth + lambda_nonCg_unmeth) * 100, 3)
    
    # Create tidy output
    cols = ["Sample_Name", "Methylated_All", "Unmethylated_All", "Percent_Methylated_All", 
            "Methylated_CpG", "Unmethylated_CpG", "Percent_Methylated_CpG", 
            "Methylated_nonCpG", "Unmethylated_nonCpG", "Percent_Methylated_nonCpG",
            "Lambda_Methylated_All", "Lambda_Unmethylated_All", "Lambda_Percent_Methylated_All", 
            "Lambda_Methylated_CpG", "Lambda_Unmethylated_CpG", "Lambda_Percent_Methylated_CpG", 
            "Lambda_Methylated_nonCpG", "Lambda_Unmethylated_nonCpG", "Lambda_Percent_Methylated_nonCpG",
            "z", "Z", "x", "X", "H", "h", "H", "u", "U", 
            "z_lambda", "Z_lambda", "x_lambda", "X_lambda", "H_lambda", "h_lambda", "H_lambda", "u_lambda", "U_lambda"]
    
    header = "\t".join(cols) + "\n"
    data = [str(sample_id), meth, unmeth, perc_meth,
            cg_meth, cg_unmeth, perc_cg_meth,
            nonCg_meth, nonCg_unmeth, perc_nonCg_meth,
            lambda_meth, lambda_unmeth, perc_lambda_meth,
            lambda_cg_meth, lambda_cg_unmeth, perc_lambda_cg_meth,
            lambda_nonCg_meth, lambda_nonCg_unmeth, perc_lambda_nonCg_meth,
            counts["z"], counts["Z"], counts["x"], counts["X"], counts["H"], counts["h"], counts["H"], counts["u"], counts["U"],
            lambda_counts["z"], lambda_counts["Z"], lambda_counts["x"], lambda_counts["X"], lambda_counts["H"], lambda_counts["h"], lambda_counts["H"], lambda_counts["u"], lambda_counts["U"]]
    data = "\t".join([str(x) for x in data]) + "\n"

    with open(str(outfile), "w") as out:
        out.write(header)
        out.write(data)
    /$
}