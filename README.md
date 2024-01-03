<a href="https://www.coriell.org">
    <img src="https://raw.githubusercontent.com/Apompetti-Cori/Bisulfite_Nextflow/main/multiqc_logo/Coriell_Logo.png" alt="Coriell Logo" width="500">
</a>

# Bisulfite Nextflow Pipeline

## Step 1: Clone pipeline into your project directory that contains the fastq files.
`git clone https://github.com/Apompetti-Cori/Bisulfite_Nextflow.git .`

You will now have a folder in your directory called `Bisulfite_Nextflow`.

## Step 2: Create csv file explaining the samples
- Use `sample_table_template.csv` as a guide for your sample table.
- First column should be sample containing sample ID's
- Consequent columns should be r1_L1, r1_L2, r1_L3, r1_L4, r2_L1, r2_L2, r2_L3, r2_L4
  - Example: a single end fastq should have a sample ID and a file name inhabiting the r1_L1 column
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns
  - Example: a paired end fastq should have a sample ID and file names inhabiting the r1_L1 and r2_L1 columns
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns and r2_L* columns
  - Essentially multilane files will be concatenated together into their respective reads and fed into the pipleine accordinly
   
The csv file will be passed into the pipeline using the `--sample_table` flag. Make sure to pass in the absolute path of the sample table.

## Step 3: Choose genome that the samples will be aligned to
- This step requires you to have already downloaded your reference geneome and prepared it using `bismark_genome_preparation`
- Once prepared, you can provide the folder of the genome to the `--db` flag. Make sure to pass in the absolute path of the genome folder.

## Step 4: Run pipeline
Once you have your sample table and genome added to your directory, you can run the pipline using Nextflow.

Example run: `nextflow run ./Bisulfite_Nexflow/main.nf -resume -with-report -with-trace -with-dag flowchart.mmd --db path_to_genome --multiqc_report_title Multiqc_Report_Title --sample_table path_to_sample_table`

I usually always include the -resume flag when starting a pipeline run since it will just say nothing to be resumed if a cache is not found.

## Additional flags
`--rrbs`: Passes flags to certain tools when the data is Reduced Represetation Bisulfite Sequencing (RRBS).
`--multiqc_report_title`: Title your multiqc report to be a bit more descriptive. Make sure not to include spaces in the title. 

## Dependencies:

## Tasklist:
- [ ] List dependencies
- [ ] Create Singularity/Docker Container file for the pipeline
