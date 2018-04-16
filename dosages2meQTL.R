#converts from dosage to meQTL
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
samples_w_exp <- as.character(fread("/home/angela/plosgen_revisions/FHS/samples_w_exp.txt")[1,])
samples_all <- as.character(fread("/home/wheelerlab3/Data/FHS_dosages/samples.txt")$V1)

for(i in 1:22){
  chr <- fread("zcat /home/wheelerlab3/Data/FHS_dosages/chr" %&% i %&% ".maf0.01.r20.8.dosage.txt.gz")
  colnames(chr) <- c("chr", "rsid", "bp", "A1", "A2", "MAF", samples_all)
  chr <- subset(chr, select = c("chr", "rsid", "bp", "A1", "A2", "MAF", samples_w_exp))
  
  SNP_file <- subset(chr, select = c("rsid", samples_w_exp))
  SNP_loc <- subset(chr, select = c("rsid", "chr", "bp"))
  fwrite(SNP_file, "/home/angela/plosgen_revisions/FHS/FHS_SNPs_chr" %&% i %&% ".txt", sep = "\t", col.names = T, na = "NA")
  fwrite(SNP_loc, "/home/angela/plosgen_revisions/FHS/FHS_SNPs_loc_chr" %&% i %&% ".txt", sep = "\t", col.names = T, na = "NA")
}
