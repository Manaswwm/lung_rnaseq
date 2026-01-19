#### this script will analyse the dge information passed to me by Martin for his project ####
## the goal of this script is to give primary dge outputs (heatmap and volcano plot) ##
## the input information is in the form equivalent to dge_tophits dataframe and expr_mtx matrix ##

##importing relevant packages
library(ggplot2)
library(tidyverse)
library(tmod)
library(gplots)
library(plotly)
library(pheatmap)
library(readxl)
library(stringr)
library(cowplot)
library(ggrepel)
library(biomaRt)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(tidyverse)

#declaring the biomart for humans
# hsap_mart = useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

##importing relevant functions
source("step1_makeVolcanoHeatplots.R")
source("step2_HLCAEnrichment.R")
source("step3_GeneAnnotate.R")
source("step4_PathwayAnnotate.R")

#### need to run this only once for HLCA ####

# ## reading the excel sheets that contain the cell-type gene markers from HLCA (Supp. table 4) ##
# 
# #reading in the sheet names from the excel workbook - 101 sheets
# sheets = excel_sheets("input_files/hlca/41586_2020_2922_MOESM6_ESM (1).xlsx")
# 
# #subsetting to keep only the SS2 sheets (SmartSeq 2)
# sheets = sheets[grep(x = sheets, pattern = "\\(\\SS2\\)")]
# 
# #sequentially reading each sheet and storing as a list
# hlca_data = lapply(sheets, function(X) read_excel("input_files/hlca/41586_2020_2922_MOESM6_ESM (1).xlsx", sheet = X))
# 
# #### writing a function that will go over all sheets (cell types) and list genes specific to cell types ####
# hlca_dataset_process = function(hlca_df){
# 
#   #print(hlca_df)
#   
#   #extracting the inner dataframe
#   hlca_subset = hlca_df
#   
#   #identifying the cluster type
#   clust_type = colnames(hlca_subset)[1]
#   
#   #making this into a dataframe
#   hlca_subset = as.data.frame(hlca_subset)
#   
#   #making first row as column names\
#   colnames(hlca_subset) = hlca_subset[1,]
#   
#   #removing first row
#   hlca_subset = hlca_subset[-1,]
#   
#   ## filtering for sanity checks ##
#   ## since each cluster is made up of clusters - filtering on how many cells IN and OUT 
#   ## of the cluster can express this gene, the former has to be max (>80%) and latter 
#   ## has to be min ( < 20%) -- these proportions could be made more stringent of lenient
#   
#   #making the proportions numeric
#   hlca_subset$pct_in_cluster = as.numeric(hlca_subset$pct_in_cluster)
#   hlca_subset$pct_out_cluster = as.numeric(hlca_subset$pct_out_cluster)
#   
#   #filtering for pct_in_cluster and pct_out_cluster
#   hlca_subset_filter = hlca_subset[hlca_subset$pct_in_cluster > 0.6 & hlca_subset$pct_out_cluster < 0.4,]
#   
#   #@dimension check
#   print(paste("For cluster - ", clust_type, "; I started with - ", dim(hlca_subset)[1], " and post-filtering I have - ", dim(hlca_subset_filter)[1], sep = ""))
#   
#   #checking if there is anything left to add
#   if(dim(hlca_subset_filter)[1] > 0){
#     
#     #keeping only things that I need
#     hlca_subset_filter = hlca_subset_filter[,c(1,2)]
#     
#     #adding cluster identifier
#     hlca_subset_filter$cluster = rep(clust_type)
#     
#     #return
#   return(hlca_subset_filter)
#   }else{return(NULL)}
# }
# 
# #reading and filtering all the cluster tables
# hlca_subset = do.call(rbind, lapply(hlca_data, function(x){hlca_dataset_process(hlca_df = x)}))
# 
# #changing class for lfc to numeric
# hlca_subset$avg_logFC = round(as.numeric(hlca_subset$avg_logFC), 2)
# 
# #writing to a file
# write.table(file = "input_files/hlca/cluster_markers_filtered.txt", x = hlca_subset, 
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

#### reading HLCA file ####

#reading this file
hlca_subset = read.delim("input_files/hlca/cluster_markers_filtered.txt", sep = "\t", header = TRUE)

#taking cell types of interest
cell_types = unique(hlca_subset$cluster)[c(seq(1,9,by=1), 17, 19)]

#subsetting to keep only specific cell types
hlca_subset = hlca_subset[hlca_subset$cluster %in% cell_types,]

## done reading the excel sheet and identifying cluster-specific markersw ##

#declaring the file paths seperately for dge_tophits and expr_mtx
dge_tophits_path = list.files(path = "input_files/deg_files/", full.names = TRUE)

### writing a function that will automate annotating the DGE tophits
dge_automate = function(filepath){

  #reading the file
  dge_tophits = read.delim(file = filepath, header = TRUE, sep = "\t")
  
  #removing .all from the string names
  comp = strsplit(filepath, "[/.]")[[1]][4]
  
  #removing entries that are termed "NewGene" - do not have any annotations
  newgene_entries = dge_tophits$Symbol[grep(pattern = "NewGene", dge_tophits$Symbol)]
  dge_tophits = dge_tophits[!dge_tophits$Symbol %in% newgene_entries,]
  
  #removing duplicated entries - example AKAP17A, AMSTL etc - they mostly have similar lfc and p-value
  #there are 20 genes that are duplicated - only taking unique
  dge_tophits = dge_tophits[!duplicated(dge_tophits$Symbol),]
  
  #keeping only significant hits
  dge_tophits = dge_tophits[dge_tophits$Pvalue < 0.05,]
  
  ### splitting this into two ###
  
  ## first - dgeinfo ##
  dge_info = dge_tophits[,c("Symbol", "Pvalue", "log2FC")]
  
  #changing column names
  colnames(dge_info) = c("gene_symbol", "p_val", "lfc")
  
  ## second - expression matrix ##
  
  #extracting column names with FPKMs
  fpkm_cols = colnames(dge_tophits)[grep(x = colnames(dge_tophits), pattern = "_FPKM")]
  
  #extracting expression matrix using this information
  expr_mtx = dge_tophits[,fpkm_cols]
  rownames(expr_mtx) = dge_tophits$Symbol
  
  #removing FPKM from the column IDs
  colnames(expr_mtx) = str_remove(string = colnames(expr_mtx), pattern = "_FPKM")
  
  
  #setting colnames order
  col_order = c(paste(cond1,"4", sep = ""), paste(cond1,"5", sep = ""), paste(cond1,"6", sep = ""),
                paste(cond2,"4", sep = ""), paste(cond2,"5", sep = ""), paste(cond2,"6", sep = ""))

  #setting expr_mtx column names order
  expr_mtx = expr_mtx[,col_order]

  #### performing initial de genes, volcano plot and heatmap analyses ####
  step1_makeVolcanoHeatplots(dge_info = dge_info, expr_mtx = expr_mtx, comp = comp)

  #### performing enrichment analysis ####
  step2_HLCAEnrichment(dge_info = dge_info, comp = comp, hlca_subset = hlca_subset)

  #### writing a table with de genes and their annotations ####
  step3_GeneAnnotate(dge_info = dge_info, comp = comp)
  
  #### performing analysis that is focused on identifying top 5 KEGG and BP terms and their terms ####
  step4_PathwayAnnotate(dge_info = dge_info, comp = comp)
  
  
}
#sequentinally sending the filepaths to the function
lapply(dge_tophits_path, function(x){dge_automate(x)})
