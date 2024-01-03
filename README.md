# Bisulfite Nextflow Pipeline

## Step 1: Create csv file explaining the samples
- First column should be sample containing sample ID's
- Consequent columns should be r1_L1, r1_L2, r1_L3, r1_L4, r2_L1, r2_L2, r2_L3, r2_L4
  - Example: a single end fastq should have a sample ID and a file name inhabiting the r1_L1 column
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns
  - Example: a paired end fastq should have a sample ID and file names inhabiting the r1_L1 and r2_L1 columns
    - If it is multilane, it should have the other lanes inhabiting the other r1_L* columns and r2_L* columns
  - Essentially multilane files will be concatenated together into their respective reads and fed into the pipleine accordinly
   
The csv file will be passed into the pipeline using the `--sample_table` flag. Make sure to pass in the absolute path of the sample table.
