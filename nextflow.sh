export NXF_WORK=$(pwd)/.work/

nextflow -log $(pwd)/.logs/nextflow.log \
run $(find $(pwd) -maxdepth 2 -type f -name "test.nf") \
-resume \
-with-report $(pwd)/.logs/run-report.html \
-with-trace $(pwd)/.logs/trace.txt \
--genome "mm9lambda" \
--multiqc_report_title "Test Title" \
--sample_table $(find $(pwd) -type f -name "sample_table.csv") \
--input_type "fastq"