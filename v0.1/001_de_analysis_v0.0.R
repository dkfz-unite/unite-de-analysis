#### Declare libraries ####
library(DESeq2)           # differential expression analysis

#### Inputs
args = commandArgs(trailingOnly = T)

#### Count matrix
count_raw = read.table(file = args[1], header = T, sep = "\t", check.names = F)
rownames(count_raw) = count_raw$gene_id
gene_tab = count_raw[,1:2]
row.names(gene_tab) = gene_tab[,1]
count_raw = count_raw[,c(-1,-2)]

#### Metadata
metadata = read.table(file = args[2], header = T, sep = "\t", check.names = F)
rownames(metadata) = metadata[,1]
metadata[,2] = as.factor(metadata[,2])

#### Count matrix re-arrangement
count_raw = count_raw[,rownames(metadata)]

#structure(metadata)

#### DESeq2
# Object initialization
dds = DESeqDataSetFromMatrix(countData = count_raw,
                             colData = metadata,
                             design = ~ condition)

# Pre-filtering: keep rows which more than 50% columns that have counts >= 10
keep = rowSums(counts(dds) >= 10) >= (ncol(count_raw)/2)
dds = dds[keep,]

# DE analysis
dds = DESeq(dds)
res = results(dds)

# Normalized matrix
count_norm = counts(dds, normalized=T)

res_tab = as.data.frame(na.omit(res))
row.names(res_tab) = res_tab$ID
res_tab$gene_name = gene_tab[row.names(res_tab),2]

write.table(x = data.frame("ID" = rownames(res_tab), res_tab), file = args[3], sep = "\t", quote = F, row.names = F)
