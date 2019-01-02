#extract snp info (pos, ref, alt) for joining to eqtl data
import itertools
import gzip

eqtl_anno = open("eqtl_anno.txt", "w")
eqtl_anno.write("snps,pos_snps,ref,alt\n")
#with gzip.open("test.vcf.gz") as f:
with gzip.open("All_20180423.vcf.gz") as f:
    for line in itertools.islice(f, 56, None): #skip first 56 lines (notes)
        snp_info = line.split("\t")
        (CHROM, POS, ID, REF, ALT) = snp_info[0:5] 
        if len(REF) == 1 and len(ALT) == 1: #keep only SNPs, no indels
            eqtl_anno.write(ID + "," + POS + "," + REF + "," + ALT + "\n")
eqtl_anno.close()
