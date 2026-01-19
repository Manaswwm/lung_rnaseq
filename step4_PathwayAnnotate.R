#### using clusterprofiler to identify BP and KEGG pathway enrichment ####
step4_PathwayAnnotate = function(dge_info = dge_info, comp = comp){
  
  #arranging rankings together for upregulated and downregulated
  rankings = sign(dge_info$lfc)*(-log10(dge_info$p_val))
  names(rankings) = dge_info$gene_symbol
  rankings = sort(rankings, decreasing = TRUE)
  
  #### Hallmark genes ####
  
  #setting the collection for GSEA
  gene_set_df = msigdbr_collections(db_species = "HS")
  kegg_pathway = msigdbr(species = "Homo sapiens", collection = "H")
  kegg_genesets = kegg_pathway %>% split(x = .$gene_symbol, f = .$gs_name)
  
  #performing fgsea
  fgsea_hallmark = fgsea(pathways = kegg_genesets, stats = rankings)
  fgsea_hallmark = fgsea_hallmark[fgsea_hallmark$padj < 0.05,]
  
  #sorting
  fgsea_hallmark = fgsea_hallmark %>% arrange(desc(NES))
  
  #removing "HALLMARK_" from pathway names
  fgsea_hallmark$pathway = str_remove(fgsea_hallmark$pathway, pattern = "HALLMARK_")
  
  #converting to a dataframe
  fgsea_hallmark = as.data.frame(fgsea_hallmark)
  
  #converting the list of elements in the leadingEdge column
  fgsea_hallmark$leadingEdge = unlist(lapply(fgsea_hallmark$leadingEdge, function(x){paste(unlist(x), collapse = ",")}))
  
  #plotting hallmark sets
  fgsea_hallmark_plot = ggplot(data = fgsea_hallmark, aes(x = NES, y = reorder(pathway, NES)))+
    geom_bar(stat="identity", fill = "#c1121f")+
    labs(x = "Norm. Enrichment Score", y = "Hallmark Pathways", title = "Hallmark set", subtitle = comp)+ theme_bw()+
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10), 
          plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5, size = 14))
  
  #saving to file the output plot
  ggsave(plot = fgsea_hallmark_plot, filename = paste("fgsea_outputs/plots/", comp, "_fgsea_hallmark.png", sep = ""),
         dpi = 600, device = "png", width = 6, height = 8)
  
  #saving to file the output table
  write.table(x = fgsea_hallmark, file = paste("fgsea_outputs/tables/", comp, "_fgsea_hallmark_table.txt", sep = ""),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  #### Biological processes ####
  
  #setting the collection for GSEA
  #gene_set_df = msigdbr_collections(db_species = "HS")
  bp_pathway = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  bp_genesets = bp_pathway %>% split(x = .$gene_symbol, f = .$gs_name)

  #performing fgsea
  fgsea_bp = fgsea(pathways = bp_genesets, stats = rankings)
  fgsea_bp = fgsea_bp[fgsea_bp$padj < 0.05,]
  
  #sorting
  fgsea_bp = fgsea_bp %>% arrange(desc(NES))
  
  #removing "HALLMARK_" from pathway names
  fgsea_bp$pathway = str_remove(fgsea_bp$pathway, pattern = "GOBP_")
  
  #converting to a dataframe
  fgsea_bp = as.data.frame(fgsea_bp)
  
  #taking only top 10
  fgsea_bp = fgsea_bp[1:10,]
  
  #converting the list of elements in the leadingEdge column
  fgsea_bp$leadingEdge = unlist(lapply(fgsea_bp$leadingEdge, function(x){paste(unlist(x), collapse = ",")}))
  
  #plotting hallmark sets
  fgsea_bp_plot = ggplot(data = fgsea_bp, aes(x = NES, y = reorder(pathway, NES)))+
    geom_bar(stat="identity", fill = "#c1121f")+
    labs(x = "Norm. Enrichment Score", y = "Processes", title = "Biological Processes", subtitle = comp)+ theme_bw()+
    theme(axis.title = element_text(size = 14), axis.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5, size = 14))
  
  #saving to file the output plot
  ggsave(plot = fgsea_bp_plot, filename = paste("fgsea_outputs/plots/", comp, "_fgsea_bp.png", sep = ""),
         dpi = 600, device = "png", width = 6, height = 8)
  
  #saving to file the output table
  write.table(x = fgsea_bp, file = paste("fgsea_outputs/tables/", comp, "_fgsea_bp_table.txt", sep = ""),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

}

## seperating the upregulated and downregulated dge ##

# #upregulated genes
# upreg_info = dge_info[dge_info$lfc > 0,]
# 
# #downregulated genes
# downreg_info = dge_info[dge_info$lfc < 0,]
# 
# ## converting gene symbol to entrez id for kegg analysis ##
# 
# #upregulated genes
# entrez_upreg = getBM(attributes = c("entrezgene_id", "external_gene_name"), filters = "external_gene_name",
#                             values = upreg_info$gene_symbol, mart = hsap_mart)
# 
# #adding entrez id
# entrez_upreg_add = merge(upreg_info, entrez_upreg, by.x = "gene_symbol", by.y = "external_gene_name")
# 
# #sorting (ranking) on lfc
# entrez_upreg_add = entrez_upreg_add %>% arrange(desc(lfc))
# 
# #downregulated genes
# entrez_downreg = getBM(attributes = c("entrezgene_id", "external_gene_name"), filters = "external_gene_name",
#                      values = downreg_info$gene_symbol, mart = hsap_mart)
# 
# #adding entrez id
# entrez_downreg_add = merge(downreg_info, entrez_downreg, by.x = "gene_symbol", by.y = "external_gene_name")
# 
# #sorting (ranking) on lfc
# entrez_downreg_add = entrez_downreg_add %>% arrange(lfc)
# 
# ### KEGG analysis ###
# 
# ## for up-reg genes ##
# 
# #extracting the enriched KEGG pathways 
# upreg_kegg = enrichKEGG(gene = c(entrez_upreg_add$entrezgene_id, entrez_downreg_add$entrezgene_id), organism = "hsa", keyType = "kegg")
# 
# #converting to a dataframe
# upreg_kegg = as.data.frame(upreg_kegg)
# upreg_kegg[1:5,]
# 
# 
# x = enrichGO(gene = dge_info$gene_symbol[dge_info$lfc > 0], ont = "BP", keyType = "SYMBOL", OrgDb = org.Hs.eg.db)
# 
# x = as.data.frame(x)
# 
# #for KEGG
# y = enrichKEGG(gene = unique(gene_symbol_convert$entrezgene_id), organism = "hsa", keyType = "kegg")
# 
# y = as.data.frame(y)
