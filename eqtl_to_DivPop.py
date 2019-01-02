#convert matrix eqtl format to cis-eQTL summary statistics
import pandas as pd
combined_pops = ["AFHI", "ALL"]
chrs = range(1, 23)

'''
pop = "AFHI"
chr = 22
'''

#header: snps    gene    statistic       pvalue  FDR     beta    chr     gene_name       start   end     gene_type       pos_snps        ref     alt
gene_anno = pd.read_csv("gene_annotation.txt") #gene-level data
eqtl_anno = pd.read_csv("eqtl_anno.txt") #snp-level data

for pop in combined_pops:
    DivPop = pd.DataFrame(columns=['snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta']) #to collect results from all chrs.
    for chr in chrs:
        meqtl = pd.read_csv(pop + "_Nk_10_PFs_chr" + str(chr) + "pcs_3.meqtl.cis.2018-05-12.txt", delim_whitespace = True) #snps    gene    statistic       pvalue  FDR     beta
        meqtl = meqtl.loc[meqtl['pvalue'] <= 1e-8] #sig cutoff for snp = 1e-8
        DivPop = pd.concat([DivPop, meqtl]) #add to overall pop
        print("Completed with " + pop + ", chr. " + str(chr) + ".")
    DivPop = DivPop.reset_index(drop = True) #to prevent multiple occurrances of the same index
    DivPop['gene'] = DivPop['gene'].str.replace("\.[^.]*$", "") #remove decimal and after in gene - https://stackoverflow.com/questions/19710898/regex-to-remove-everything-after-the-last-dot-in-a-file
    DivPop = DivPop.merge(gene_anno, on = "gene") #add gene data
    DivPop = DivPop.merge(eqtl_anno, on = "snps") #add snp data
    DivPop.to_csv(pop + "_cis_eqtl_summary_statistics.txt.gz", index = False, compression = "gzip", sep = "\t", na_rep = "NA")
    print("Completed with " + pop + ".")
    
