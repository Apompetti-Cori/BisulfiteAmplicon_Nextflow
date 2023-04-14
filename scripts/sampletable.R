setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sampletable <- read.csv("../sample_table_template.csv")
datapath <- "/mnt/data/research_data/2021-12-07_mouse_rrbs_peace/usftp21.novogene.com/raw_data/"
dirlist <- list.dirs(datapath,full.names = FALSE)
dirlist <- dirlist[dirlist != "" & dirlist != "Undetermined"]
list.files(paste(datapath, dirlist,sep = ""))
sampletable[1:length(dirlist),] <- NA
sampletable$sample <- dirlist
rownames(sampletable) <- sampletable$sample
for(sample in dirlist){
  list.files(path = paste(datapath,sample,sep = "/"), pattern = "*.gz")
  fastqfiles <- list.files(path = paste(datapath,sample,sep = "/"), pattern = "*.gz")
  fastqfiles <- paste("/mnt/data/data_ap/14Apr2023_RRBS_20211207/fastq",fastqfiles,sep = "/")
  sampletable[sample,c("r1_L1","r2_L1")] <- c(fastqfiles)
}
sampletable[is.na(sampletable)] <- ""
write.csv(sampletable, file = "../sample_table_input.csv", quote = FALSE, row.names = FALSE)
