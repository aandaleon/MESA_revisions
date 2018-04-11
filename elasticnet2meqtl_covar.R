"%&%" = function(a,b) paste(a,b,sep="")
library(data.table)

for(i in c(3, 5, 10)){
  for(pop in c("AFA", "CAU", "HIS")){
    elasticnet2meqtl_covar <- fread("/home/angela/plosgen_revisions/MESA/covariates/" %&% pop %&% "pcs" %&% i %&%"cov.txt")
    elasticnet2meqtl_covar$FID <- NULL
    elasticnet2meqtl_covar$pop <- NULL
    elasticnet2meqtl_covar$IID <- as.character(elasticnet2meqtl_covar$IID)
    elasticnet2meqtl_covar <- transpose(elasticnet2meqtl_covar)
    fwrite(elasticnet2meqtl_covar, "/home/angela/plosgen_revisions/MESA/covariates/" %&% pop %&% "pcs" %&% i %&%"cov.txt", col.names = F, row.names = T, sep = "\t", quote = F)
  }
}