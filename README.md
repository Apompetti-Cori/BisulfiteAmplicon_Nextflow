<img src="https://raw.githubusercontent.com/Apompetti-Cori/Bisulfite_Nextflow/main/multiqc_logo/Coriell_Logo.png" width="1000"> 

#Bisulfite Nextflow Pipeline

## Step 1: Create csv file explaining the samples
- First column should be sample containing sample ID's
- Consequent columns should be r1_L1, r1_L2, r1_L3, r1_L4, r2_L1, r2_L2, r2_L3, r2_L4
  - Example: a single end fastq should have a sample ID and a file name inhabiting the r1_L1 column
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns
  - Example: a paired end fastq should have a sample ID and file names inhabiting the r1_L1 and r2_L1 columns
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns and r2_L* columns
  - Essentially multilane files will be concatenated together into their respective reads and fed into the pipleine accordinly
   
The csv file will be passed into the pipeline using the `--sample_table` flag. Make sure to pass in the absolute path of the sample table.

## Step 2: Choose genome that the samples will be aligned to
- This step requires you to have already downloaded your reference geneome and prepared it using `bismark_genome_preparation`
- Once prepared, you can provide the folder of the genome to the `--db` flag. Make sure to pass in the absolute path of the genome folder.


## Additional flags
`--rrbs`: passes flags to certain tools when the data is Reduced Represetation Bisulfite Sequencing (RRBS).
