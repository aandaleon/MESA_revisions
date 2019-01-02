#extract gene information from gencode for joining to eqtl data
import itertools
import pandas as pd
annotate_eqtl_dups = []
with open("gencode.v28lift37.basic.annotation.gtf") as f:
#with open("gencode_test.txt") as f:
    for line in itertools.islice(f, 5, None): #skip first five lines
        gene_info = line.replace(";", "\t").split("\t")
        #each gene has different number of tags for some reason
        (chr, x1, x2, start, end) = gene_info[0:5]
        #(chr, index1, index2, start, end, index5, index6, index7, gene_id, index9, gene_type, gene_name) = gene_info[0:12]
        
        #so we work around it using something I googled - https://stackoverflow.com/questions/6039425/sequence-find-function-in-python
        gene_id = next(i for i in gene_info if i.startswith("gene_id \""))
        gene_type = next(i for i in gene_info if i.startswith(" gene_type \""))
        gene_name = next(i for i in gene_info if i.startswith(" gene_name \""))
        
        #convert data to format suitable for eqtl annotation (remove excess info)
        chr = chr.replace("chr", "")
        gene_id = gene_id.replace("gene_id \"", "").replace("\"", "").split('.')[0] #remove all but the actual id
        gene_type = gene_type.replace(" gene_type \"", "").replace("\"", "")
        gene_name = gene_name.replace(" gene_name \"", "").replace("\"", "")
        annotate_eqtl_dups.append([gene_id, chr, gene_name, start, end, gene_type])
annotate_eqtl_dups = pd.DataFrame(annotate_eqtl_dups) #df easier to work w/ than list of list
annotate_eqtl_dups.columns = ["gene", "chr", "gene_name", "start", "end", "gene_type"]

#keep only full gene start and end
gene_list = list(set(annotate_eqtl_dups.iloc[:,2].tolist())) #force unique
annotate_eqtl = []
for gene in gene_list:
    gene_df = annotate_eqtl_dups.loc[annotate_eqtl_dups["gene_name"] == gene]
    start = gene_df['start'].min()
    end = gene_df['end'].max()
    gene_ensg = gene_df.iloc[0]['gene']
    chr = gene_df.iloc[0]['chr']
    gene_type = gene_df.iloc[0]['gene_type']
    annotate_eqtl.append([gene_ensg, chr, gene, start, end, gene_type])
annotate_eqtl = pd.DataFrame(annotate_eqtl)
annotate_eqtl.columns = ["gene", "chr", "gene_name", "start", "end", "gene_type"]
annotate_eqtl.to_csv("gene_annotation.txt", index = False)

