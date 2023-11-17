#### Libraries ####
library(stringr)                                # String manipulation
library(clusterProfiler)                        # GSEA
library(enrichplot)                             # GSEA result visualization
library(ggplot2)                                # GSEA result visualization
library("org.Hs.eg.db", character.only = TRUE)  # Annotation for GSEA

#### Inputs ####
#1: DESeq2 result table (all)
#2: Output directory

args = commandArgs(trailingOnly = T)
res_file  = args[1] #1
out_dir   = args[2] #2

#### Main ####
res_tab = read.table(file = res_file, header = T, sep = "\t")
rownames(res_tab) = str_remove(string = res_tab[,1], pattern = '\\.[0-9]+') # Remove version ID
res_tab = res_tab[,-1]

# GSEA
gsea_log2fc = res_tab$log2FoldChange
names(gsea_log2fc) = rownames(res_tab)
gsea_log2fc = sort(gsea_log2fc, decreasing = T)

gsea_analysis = gseGO(
  geneList = gsea_log2fc,
  ont = "ALL",
  keyType = "ENSEMBL",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  verbose = T,
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "BH",
  nPermSimple = 1000000,
  eps = 0
)

gsea_dotplot = dotplot(gsea_analysis, showCategory=10, title="Enriched GO Terms", split=".sign") + facet_grid(.~.sign)
gsea_treeplot = treeplot(pairwise_termsim(gsea_analysis), fontsize = 3, offset = 16)

ggsave(filename = paste(out_dir, "gsea_dotplot.pdf", sep = "/"), plot = gsea_dotplot, device = "pdf")
ggsave(filename = paste(out_dir, "gsea_treeplot.pdf", sep = "/"), plot = gsea_treeplot, device = "pdf", width = 10)

write.table(x = gsea_analysis@result, file = paste(out_dir, "gsea_result.tsv", sep = "/"), sep = "\t", quote = F, col.names = NA)
# KEGG
entrez_ids = bitr(geneID = rownames(res_tab), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
res_tab_kegg = res_tab[entrez_ids$ENSEMBL,]
res_tab_kegg$entrezid = entrez_ids$ENTREZID

kegg_log2fc = res_tab_kegg$log2FoldChange
names(kegg_log2fc) = res_tab_kegg$entrezid
kegg_log2fc = sort(kegg_log2fc, decreasing = T)

kegg_analysis = gseKEGG(
  geneList = kegg_log2fc,
  organism = "hsa", # Homo sapiens
  nPermSimple = 1000000,
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = "ncbi-geneid"
)

kegg_dotplot = dotplot(kegg_analysis, showCategory=10, title="Enriched KEGG Pathways", split=".sign") + facet_grid(.~.sign)
kegg_emap = emapplot(pairwise_termsim(kegg_analysis))

ggsave(filename = paste(out_dir, "kegg_dotplot.pdf", sep = "/"), plot = kegg_dotplot, device = "pdf")
ggsave(filename = paste(out_dir, "kegg_emap.pdf", sep = "/"), plot = kegg_emap, device = "pdf")

write.table(x = kegg_analysis@result, file = paste(out_dir, "kegg_result.tsv", sep = "/"), sep = "\t", quote = F, col.names = NA)
