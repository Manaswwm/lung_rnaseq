#### writing a function that will write to a table the gene ontology ####
step3_GeneAnnotate = function(dge_info, comp){
  
  #sorting according to the log fold change (absolute)
  dge_info_sorted = dge_info[order(abs(dge_info$lfc), decreasing = TRUE),]
  
  #extracting the gene annotation
  gene_annotation = getBM(attributes = c("external_gene_name", "goslim_goa_description"), 
        filters = c("external_gene_name"), values = dge_info_sorted$gene_symbol, mart = hsap_mart, useCache = FALSE)
  
  #collapsing multiple gene annotations per gene into one
  gene_annotation_collapse = as.data.frame(gene_annotation %>% group_by(external_gene_name) %>% summarise(Annotation = paste(goslim_goa_description, collapse = " | ")))
  
  #changing column names for gene_symbol
  colnames(gene_annotation_collapse)[1] = "gene_symbol"
  
  #merging the gene annotation with log fold changes
  dge_info_sorted = merge(dge_info_sorted, gene_annotation_collapse, by = "gene_symbol")
  
  #resorting
  dge_info_sorted = dge_info_sorted[order(abs(dge_info_sorted$lfc), decreasing = TRUE),]

  #writing to a file
  write.table(file = paste("annotation_outputs/annotation_", comp, ".txt", sep = ""),
              x = dge_info_sorted, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  }
