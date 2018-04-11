#compares the beta values for the old mEQTL (downloaded imputation) vs new meQTL (Michigan imputation)
"%&%" = function(a,b) paste(a,b,sep="")
library(data.table)
library(ggplot2)
args <- commandArgs(trailingOnly=T)
pop <- args[1] #ex. MXL, YRI
Nk <- args[2] #number of PEER factors
pcs <- args[3]

old.meQTL <- "/home/lauren/files_for_revisions_plosgen/old_meqtl_results/"
new.meQTL <- "/home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/"
out.dir <- "/home/angela/plosgen_revisions/compare_meQTL_MESA/"

i <- 1
old <- fread("zcat " %&% old.meQTL %&% pop %&% "_Nk_" %&% Nk %&% "_PFs_chr" %&% i %&% ".meqtl.cis.2017-07-05.txt.gz")
new <- fread("zcat " %&% new.meQTL %&% pop %&% "_Nk_" %&% Nk %&% "_PFs_chr" %&% i %&% "pcs_" %&% pcs %&% ".meqtl.cis.2018-04-04.txt.gz")
old_v_new <- merge(old, new, by = c("snps", "gene"))
colnames(old_v_new) <- c("snps", "gene", "old_t", "old_p", "old_FDR", "old_beta", "new_t", "new_p", "new_FDR", "new_beta")

for(i in c(2:22)){
  #all old meQTL have 10 PCs (lines 255-256)
  old_add <- fread("zcat " %&% old.meQTL %&% pop %&% "_Nk_" %&% Nk %&% "_PFs_chr" %&% i %&% ".meqtl.cis.2017-07-05.txt.gz")
  new_add <- fread("zcat " %&% new.meQTL %&% pop %&% "_Nk_" %&% Nk %&% "_PFs_chr" %&% i %&% "pcs_" %&% pcs %&% ".meqtl.cis.2018-04-04.txt.gz")
  old_v_new_add <- merge(old, new, by = c("snps", "gene"))
  colnames(old_v_new_add) <- c("snps", "gene", "old_t", "old_p", "old_FDR", "old_beta", "new_t", "new_p", "new_FDR", "new_beta")
  old_v_new <- rbind(old_v_new, old_v_new_add)
  cat('Chromosome ' %&% i %&% ' has been added to the model successfully.')
}

#beta
r2_beta <- round((cor.test(old_v_new$old_beta, old_v_new$new_beta)$estimate)^2, 3)
png(filename = out.dir %&% "betas_" %&% pop %&% "_" %&% Nk %&% "_PFs_" %&% pcs %&% "_PCs.png", width = 550, height = 500)
ggplot() + 
  geom_point(data = old_v_new, aes(x = old_beta, y = new_beta)) + 
  labs(title = "Betas: " %&% pop %&% ", " %&% Nk %&% " PFs, " %&% pcs %&% " PCs", subtitle = "r2 = " %&% r2_beta, x = "old meQTL beta", y = "new meQTL beta") +
  coord_fixed() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

#t
r2_t <- round((cor.test(old_v_new$old_t, old_v_new$new_t)$estimate)^2, 3)
png(filename = out.dir %&% "t_" %&% pop %&% "_" %&% Nk %&% "_PFs_" %&% pcs %&% "_PCs.png", width = 550, height = 500)
ggplot() + 
  geom_point(data = old_v_new, aes(x = old_t, y = new_t)) + 
  labs(title = "t: " %&% pop %&% ", " %&% Nk %&% " PFs, " %&% pcs %&% " PCs", subtitle = "r2 = " %&% r2_t, x = "old meQTL t", y = "new meQTL t") +
  coord_fixed() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

#-log10(p)
r2_log10p <- round((cor.test(-log10(old_v_new$old_p), -log10(old_v_new$new_p))$estimate)^2, 3)
png(filename = out.dir %&% "-log10p_" %&% pop %&% "_" %&% Nk %&% "_PFs_" %&% pcs %&% "_PCs.png", width = 550, height = 500)
ggplot() + 
  geom_point(data = old_v_new, aes(x = -log10(old_p), y = -log10(new_p))) + 
  labs(title = "-log10(p): " %&% pop %&% ", " %&% Nk %&% " PFs, " %&% pcs %&% " PCs", subtitle = "r2 = " %&% r2_log10p, x = "old meQTL -log10(p)", y = "new meQTL -log10(p)") +
  coord_fixed() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
dev.off()

cat('Analysis done in: ' %&% pop %&% ", " %&% Nk %&% " PFs, " %&% pcs %&% " PCs\n")

