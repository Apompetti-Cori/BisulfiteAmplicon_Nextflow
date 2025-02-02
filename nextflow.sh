#!/usr/bin/env bash

export NXF_WORK=$(pwd)/.work/

nextflow -log $(pwd)/.logs/nextflow.log \
run $(find $(pwd) -maxdepth 2 -type f -name "main.nf") \
-resume \
-with-report $(pwd)/.logs/run-report.html \
-with-trace $(pwd)/.logs/trace.txt \
-with-dag $(pwd)/.logs/flowchart.mmd \
--genome "mm9lambda" \
--multiqc_report_title "Test Title" \
--sample_table $(find $(pwd) -type f -name "sample_table.csv") \
--input_type "fastq"