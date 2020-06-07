##arguments
compareCluster_autoCena= function(cluster_data=cluster_information){



##args end

cluster_data = cluster_data %>% dplyr::filter(cluster_included=="yes")
cluster_available = cluster_data %>% pull(clusters)

biggest_cluster = cluster_data %>% dplyr::filter(gene_no==max(cluster_data$gene_no))


gene_universe = row.names(ds)
normdata= t(ds) %>% as.data.frame()
normdata$merged = purrr::pmap(corresp_info[intersect(voi,colnames(info_dataset))], paste, sep="-") %>% unlist()



#setwd(paste0(originalwd,cutoff_wd))
#
# plotPath = paste0(originalwd, cutoff_wd, "/clusterProfiler")
# dir.create(plotPath, showWarnings = FALSE)

orga_type = tolower(organism)

#mouse code
if(orga_type == "mouse"){
  universe_mouse_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = universe, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  universe_mouse_human <- universe_mouse_human[,2]
  universe_Entrez_mouse = clusterProfiler:: bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  universe_Entrez_mouse = unlist(universe_Entrez_mouse[2],use.names = FALSE)
  universe_Entrez_mouse_human = clusterProfiler::bitr(universe_mouse_human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  universe_Entrez_mouse_human = unlist(universe_Entrez_mouse_human[2],use.names = FALSE)
}else{
#humans @shobhit choose which library function to use
universe_Entrez = clusterProfiler:: bitr(gene_universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = F)
#library('org.Hs.eg.db')
}
c1_hallmark_genes = gmtfile_hallmarks

#apply function to cluster

list_of_entrez <- list()

#rn =2

for(rn in 1:nrow(cluster_data)){
#cluster_meta_annot_excel =   function(rn, clusterdf) {
  dv= cluster_data[rn, ]
  print(paste0("Prediction for ", dv$color," containing ", dv$gene_no, " genes"))
  genes_dv = stri_split_regex(dv$gene_n, pattern=",") %>% unlist()

  ##mouse code
  if(orga_type == "mouse"){
    universe_orignal <- universe
    genes_cluster = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unlist(genes_dv), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    cluster_genes <- genes_cluster[,2]
    list_of_genes_mouse_human <- list(cluster_genes)
    entrez_de = clusterProfiler::bitr(unlist(genes_dv), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    module_entrez_mouse <- unlist(entrez_de[2],use.names = FALSE)
    entrez_de = clusterProfiler::bitr(unlist(list_of_genes_mouse_human), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    module_entrez_mouse_human <- unlist(entrez_de[2],use.names = FALSE)

    list_of_entrez[[dv$color]] <- module_entrez_mouse


  }else{

  #humans
  #@thomas shobhit some symbols do not have entrezids, maybe they have aliases  e.g. C6orf48
  #select entrez ids for the genes
  module_entrez= universe_Entrez %>%
    #drop_na(ENTREZID) %>%
    filter(SYMBOL%in%genes_dv) %>%
    pull(ENTREZID)

  list_of_entrez[[dv$color]] <- module_entrez


  }

  }



ck <- compareCluster(geneCluster = list_of_entrez, fun = "enrichGO",
                     OrgDb='org.Hs.eg.db',
                     pvalueCutoff  = 0.05,
                     pAdjustMethod = "none",
                     ont = "BP",
                     univers = universe_Entrez$ENTREZID,
                     qvalueCutoff = 1)

return(ck)

# dotplot(ck,showCategory=numTerms) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=15),
#                                     axis.text.x  = element_text(angle=90, vjust=0.5, size=10),axis.text.y  = element_text(size=10))


}
