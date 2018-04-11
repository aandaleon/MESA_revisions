"%&%" = function(a,b) paste(a,b,sep="")
library(data.table)
library(dplyr)

pop <- "YRI"

for(i in 1:22){
  chr <- fread("/home/angela/plosgen_revisions/HapMap/" %&% pop %&% "/" %&% pop %&% "_" %&% i %&%".dosage.txt")
  samples <- fread("/home/angela/plosgen_revisions/HapMap/" %&% pop %&% "/samples_" %&% pop %&% ".txt", header = F)
  colnames(chr) <- c("chr", "rsid", "pos", "A1", "A2", "MAF", samples$V1)
  samples <- samples$V1
  SNP_file <- select(chr, rsid, samples)
  snps_location <- select(chr, rsid, chr, pos)
  snps_location$chr <- paste("chr", snps_location$chr, sep = "")
  fwrite(SNP_file, "/home/angela/plosgen_revisions/HapMap/" %&% pop %&% "/SNP/chr" %&% i %&% ".txt", col.names = T, sep = "\t")
  fwrite(snps_location, "/home/angela/plosgen_revisions/HapMap/" %&% pop %&% "/snps_location/chr" %&% i %&% ".txt", col.names = T, sep = "\t")
}
