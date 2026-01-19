#### performing customized enrichment with the HLCA dataset ####
step2_HLCAEnrichment = function(dge_info, comp, hlca_subset){
  
  #sorting the dge_info df to annotate up/downregulation
  dge_info = dge_info %>% mutate(regulation = ifelse(lfc > 0, "up", ifelse(lfc < 0, "down", "neutral"))) %>%
    filter(regulation != "neutral")
    
  #identifying all the genes represented in the clusters
  cluster_genes_all = unique(hlca_subset$Gene)
  
  #### writing a function to perform cluster enrichment
  compute_enrichment = function(de_genes, gene_cluster, cluster_genes_all, type) {
      
    #subsetting to keep the cell type and the gene ID - this is done with this reference
    #https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
    term2gene = hlca_subset[,c(3,1)]
    
    #renaming the cell type and gene ID columns
    colnames(term2gene) = c("cell_type", "gene_id")
    
    #performing enrichment using enricher function directly
    enrich_res = enricher(gene = de_genes, TERM2GENE = term2gene)
    
    #plotting the output as dotplot
    enrichplot = barplot(enrich_res)
    
    enrichplot = enrichplot + ggtitle(label = type, subtitle = comp)+xlab("Gene Count")+
      theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16, angle = 30),
            plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 14, hjust = 0.5))
    
    #returing the mega tidy element that now contains the enrichment analysis
    return(enrichplot)
  }
  
  #splitting the up- and downregulated genes and analysing them seperately
  up_genes = dge_info %>% filter(regulation == "up") %>% pull(gene_symbol)
  down_genes = dge_info %>% filter(regulation == "down") %>% pull(gene_symbol)
  
  #computing the enrichment seperately
  enrich_up = compute_enrichment(de_genes = up_genes, gene_cluster = hlca_subset, 
                                 cluster_genes_all = cluster_genes_all, type = "Upregulated genes")
  enrich_down = compute_enrichment(de_genes = down_genes, gene_cluster = hlca_subset, 
                                   cluster_genes_all = cluster_genes_all, type = "Downregulated genes")
  
  #arranging
  enrich_plot = ggarrange(enrich_up, enrich_down)
  
  #saving with ggsave
  ggsave(filename = paste("output_plots/enrichment_plot_enricher/enrich_plot_", comp, ".jpeg", sep = ""), plot = enrich_plot,
         dpi = 300, height = 8, width = 14, device = "jpeg", bg = "white")
  
  #@print check
  print(paste("Gone over for enrichment analysis - ", comp, sep = ""))

}



#### un-used homemade code ####  

# #mega tidy element
# enrichment = gene_cluster %>%
#   group_by(cluster) %>%
#   summarise(
#     in_cluster = sum(Gene %in% de_genes),
#     total_in_cluster = n(),
#     .groups = "drop"
#   ) 
# 
# enrichment2 = enrichment %>%
#   mutate(
#     not_in_cluster = length(de_genes) - in_cluster,
#     total_not_in_cluster = length(cluster_genes_all) - total_in_cluster
#   ) 
# 
# enrichment3 = enrichment2 %>%
#   rowwise() %>%
#   mutate(
#     p_value = fisher.test(
#       matrix(c(in_cluster,
#                total_in_cluster - in_cluster,
#                not_in_cluster,
#                total_not_in_cluster - not_in_cluster),
#              nrow = 2)
#     )$p.value
#   )
# 
# enrichment3 = enrichment2 %>%
#   rowwise() %>%
#   mutate(
#       p_value = fisher.test(matrix(c(enrichment2$in_cluster,
#                enrichment2$total_in_cluster - enrichment2$in_cluster,
#              enrichment2$not_in_cluster,
#              abs(enrichment2$total_not_in_cluster - enrichment2$not_in_cluster)),
#              nrow = 2))$p.value)
# 
# enrichment4= enrichment3 %>%
#   ungroup() %>%
#   mutate(
#     gene_ratio = in_cluster / total_in_cluster
#   ) %>%
#   arrange(p_value)




# #extracting only the top 10 terms
# top_up = enrich_up %>% filter(p_value < 0.05) %>% arrange(desc(gene_ratio))%>% arrange(desc(gene_ratio))
# top_down = enrich_down %>% filter(p_value < 0.05) %>% arrange(desc(gene_ratio)) %>% arrange(desc(gene_ratio))
# 
# ####writing a function to plot the enrichment analysis
# plot_enrichment <- function(df, title) {
#   
#   #writing a ggplot code
#   plot = ggplot(df, aes(x = fct_reorder(cluster, gene_ratio),
#                  y = -log10(p_value),
#                  size = gene_ratio,
#                  color = p_value)) +
#     geom_point() +
#     coord_flip() +
#     scale_color_gradient(low = "red", high = "blue", name = "p-value") +
#     labs(
#       x = "Cluster",
#       y = "-log10(p-value)",
#       title = paste(title, " for - ", comp, sep = "")
#     ) +
#     theme_minimal(base_size = 14)+
#     theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text = element_text(size = 14), 
#           axis.title = element_text(size = 16), legend.title = element_text(size = 16),
#           legend.text = element_text(size = 14))
#   
#   #returning
#   return(plot)
# }
# 
# #plotting the up and downreg enrichment seperately
# plot_up = plot_enrichment(top_up, "Upreg genes")
# plot_down = plot_enrichment(top_down, "Downreg genes")