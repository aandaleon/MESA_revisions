#make PLINK-input dosages from PrediXcan and then make PC files using smartpca 
"%&%" = function(a,b) paste(a,b,sep="")
library(data.table)
in.dir <- "/home/wheelerlab3/Data/FHS_dosages/"
out.dir <- "/home/angela/plosgen_revisions/FHS/PCA/"
pop <- "FHS"

fam <- fread(in.dir %&% "samples.txt", header = F)
fam$V3 <- 0
fam$V4 <- 0
fam$V5 <- 0
fam$V6 <- 1 #smartpca requires a phenotype
fwrite(fam, out.dir %&% pop %&% ".fam", sep = "\t", col.names = F)

#concatenated all dosage files b/c plink2 only takes 1 at a time and that's annoying
system("/home/angela/plosgen_revisions/plink2 --fam " %&% out.dir %&% pop %&% ".fam --import-dosage " %&% out.dir %&% pop %&%".dosage.txt.gz skip0=1 skip1=1 chr-col-num=1 pos-col-num=3 noheader format=1 --make-bed --out " %&% out.dir %&% pop)
system("plink --bfile " %&% out.dir %&% pop %&% " --recode --out " %&% out.dir %&% pop)

#LD-prune
system("/usr/local/bin/plink --bfile " %&% out.dir %&% pop %&% " --indep-pairwise 50 5 0.2 --out " %&% out.dir %&% "pruned")

#Extract LD-pruned and keep only SNPs
system("/usr/local/bin/plink --bfile " %&% out.dir %&% pop %&% " --extract " %&% out.dir %&% "pruned.prune.in --snps-only --make-bed --out ld_pruned")

#Make map and ped
system("/usr/local/bin/plink --bfile " %&% out.dir %&% "ld_pruned --recode --out " %&% out.dir %&% "ld_pruned")

#Run smartpca
system("/usr/bin/smartpca -p " %&% out.dir %&% pop %&% "_smartpca.par")

#Make PCs into HapMap format
PCs <- fread(out.dir %&% "ld_pruned.evec", sep = " ")
PCs$V1 <- gsub(":.*", "", PCs$V1) #makes first column just FID
PCs$V12 <- NULL
colnames(PCs) <- c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
PC_names <- transpose(as.data.frame(colnames(PCs)))
colnames(PC_names) <- colnames(PCs)
PCs <- transpose(rbind(PC_names, PCs))
fwrite(PCs, out.dir %&% pop %&% "_10_PCs.txt", col.names = F, na = "NA", sep = "\t")
PCs_5 <- PCs[1:6,]
fwrite(PCs_5, out.dir %&% pop %&% "_5_PCs.txt", col.names = F, na = "NA", sep = "\t")
PCs_3 <- PCs[1:4,]
fwrite(PCs_5, out.dir %&% pop %&% "_3_PCs.txt", col.names = F, na = "NA", sep = "\t")


