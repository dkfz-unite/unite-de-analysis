#### Declare libraries ####
library(DESeq2)       # DE analysis
library(ggplot2)      # Plotting
library(plotly)       # Plotting
library(htmlwidgets)  # Export plotly to html format

#### Inputs ####
#1: raw count matrix: in tsv format, must have "gene_id" & "gene_name" columns as the 1st two columns
#2: metadata table: "col_id" column is the column names of raw count matrix
#3: name of the "in-comparison" column (in metadata table): conditions, treatments, cell lines, etc.
#4: output directory
#5: mode: "simple" or "pairwise"
#6: "control" sample type: the name of control sample type in #3 column. Only applicable for "simple" mode

args = commandArgs(trailingOnly = T)
count_raw_file  = args[1] #1
metadata_file   = args[2] #2
compar_col      = args[3] #3
output_dir      = args[4] #4
run_mode        = args[5] #5
control_name    = args[6] #6

#### Parameters check ####

#### Pre-processing ####
# Count matrix
count_raw = read.table(file = count_raw_file, header = T, sep = "\t", check.names = F)
rownames(count_raw) = count_raw$gene_id
gene_ref = count_raw[,1:2]
count_raw = count_raw[,c(-1,-2)]

# Metadata
metadata = read.table(file = metadata_file, header = T, sep = "\t", check.names = F)
rownames(metadata) = metadata[,"col_id"]
metadata[,compar_col] = as.factor(metadata[,compar_col])

# Count matrix re-arrangement
count_raw = count_raw[,rownames(metadata)]

#### DESeq2 ####
# Object initialization
design_str = paste0('~ ',compar_col)
dds = DESeqDataSetFromMatrix(countData = count_raw,
                             colData = metadata,
                             design = formula(design_str))

# Pre-filtering: keep rows which more than 50% columns that have counts >= 10
keep = rowSums(counts(dds) >= 10) >= (ncol(count_raw)/2)
dds = dds[keep,]

# DE analysis
dds = DESeq(dds)

# Normalized matrix
count_norm = counts(dds, normalized=T)

#### Volcano plot function ####
func_volcano_plot = function(res, padj, main_title, output) {
  # Pre-processing
  res = na.omit(res)
  res$consequence = ifelse(res$padj < padj,
                           ifelse(res$log2FoldChange < 0, "Down", "Up"),
                           "Not Significant")
  
  # Plotting
  temp_plot = ggplot(data = res) +
    aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name) +
    geom_point(aes(color = consequence)) +
    scale_color_manual(values = c("Up"="firebrick3", "Down"="royalblue3", "Not Significant"="grey")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    ylab("-log10(p-value)") +
    ggtitle(main_title)
  
  temp_plotly = ggplotly(temp_plot)
  saveWidget(widget = temp_plotly, file = output, selfcontained = F)
}

#### "simple" mode ####
if (run_mode == "simple") {
  # Loop over compar_col and extract all comparisons against control
  compar_col_value = metadata[,compar_col]
  for (treatment in levels(compar_col_value)) {
    if (treatment == control_name) {next}
    res = results(
      object = dds,
      contrast = c(compar_col,treatment,control_name),
      alpha = 0.05
    )
    
    # Write full result
    res_all = na.omit(as.data.frame(res)) # Remove NA rows due to DESeq2 outliers filtering
    res_all$gene_name = gene_ref[rownames(res_all),c("gene_name")]
    write.table(x = res_all,
                file = paste0(output_dir,"/deseq_",treatment,"_",control_name,"_all.tsv"),
                sep = "\t",
                col.names = NA
    )
    
    # Write up-regulated hits
    res_up = res_all[res_all$padj < 0.05 & res_all$log2FoldChange > 0,]
    write.table(x = res_up,
                file = paste0(output_dir,"/deseq_",treatment,"_",control_name,"_up.tsv"),
                sep = "\t",
                col.names = NA
    )
    
    # Write down-regulated hits
    res_down = res_all[res_all$padj < 0.05 & res_all$log2FoldChange < 0,]
    write.table(x = res_down,
                file = paste0(output_dir,"/deseq_",treatment,"_",control_name,"_down.tsv"),
                sep = "\t",
                col.names = NA
    )
    
    # Volcano plot
    volcano_output = paste0(output_dir,"/deseq_",treatment,"_",control_name,"_volcano.html")
    volcano_title = paste0(toupper(treatment)," vs. ",toupper(control_name))
    func_volcano_plot(res = res_all, padj = 0.05, main_title = volcano_title, output = volcano_output)
  }
} else {
  #### "pairwise" mode ####  
  compar_col_value = metadata[,compar_col]
  for (i in 1:(length(levels(compar_col_value))-1)) {
    cond_i = levels(compar_col_value)[i]
    for (j in (i+1):length(levels(compar_col_value))) {
      cond_j = levels(compar_col_value)[j]
      
      # Extract result
      contrast = c(cond_i, cond_j)
      res = results(
        object = dds,
        contrast = c(compar_col, contrast),
        alpha = 0.05
      )
      
      # Write full result
      res_all = na.omit(as.data.frame(res)) # Remove NA rows due to DESeq2 outliers filtering
      res_all$gene_name = gene_ref[rownames(res_all),c("gene_name")]
      write.table(x = res_all,
                  file = paste0(output_dir,"/deseq_",cond_i,"_",cond_j,"_all.tsv"),
                  sep = "\t",
                  col.names = NA
      )
      
      # Write up-regulated hits
      res_up = res_all[res_all$padj < 0.05 & res_all$log2FoldChange > 0,]
      write.table(x = res_up,
                  file = paste0(output_dir,"/deseq_",cond_i,"_",cond_j,"_up.tsv"),
                  sep = "\t",
                  col.names = NA
      )
      
      # Write down-regulated hits
      res_down = res_all[res_all$padj < 0.05 & res_all$log2FoldChange < 0,]
      write.table(x = res_down,
                  file = paste0(output_dir,"/deseq_",cond_i,"_",cond_j,"_down.tsv"),
                  sep = "\t",
                  col.names = NA
      )
      
      # Volcano plot
      volcano_output = paste0(output_dir,"/deseq_",cond_i,"_",cond_j,"_volcano.html")
      volcano_title = paste0(toupper(cond_i)," vs. ",toupper(cond_j))
      func_volcano_plot(res = res_all, padj = 0.05, main_title = volcano_title, output = volcano_output)
    }
  }
}
