setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sampletable <- read.csv("../sample_table_input.csv")
row.names(sampletable) <- sampletable$sample
datapath <- "/mnt/data/research_data/2021-03-29_rrbs61-65_peace_hg_himani_mm/usftp21.novogene.com/raw_data/"
for(sample in sampletable$sample){
  print(list.files(path = paste(datapath,sample,sep = ""), pattern = "*.gz"))
  fastqfiles <- list.files(path = paste(datapath,sample,sep = ""), pattern = "*.gz")
  fastqfiles <- paste("/home/apompetti/data/13Apr2023_RRBS_20210329/fastq/",fastqfiles,sep = "")
  sampletable[sample,c("r1_L1","r2_L1")] <- fastqfiles
}
