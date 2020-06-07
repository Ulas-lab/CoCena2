##arguments
clusterprofiler_autoCena= function(cluster_data=cluster_information,
                                   cutoff_wd=cutoff_wd,
                                   originalwd=working_directory,
                                   chosen_cutoff=optimal_cutoff,
                                   group=voi){



##args end

cluster_data = cluster_data %>% dplyr::filter(cluster_included=="yes")
cluster_available = cluster_data %>% pull(clusters)

biggest_cluster = cluster_data %>% dplyr::filter(gene_no==max(cluster_data$gene_no))

corresp_info = info_dataset[rownames(dd2)%in%rownames(info_dataset),]
gene_universe = row.names(ds)
normdata= t(ds) %>% as.data.frame()
normdata$merged = purrr::pmap(corresp_info[intersect(group,colnames(info_dataset))], paste, sep="-") %>% unlist()



#setwd(paste0(originalwd,cutoff_wd))

plotPath = paste0(originalwd, cutoff_wd, "/clusterProfiler")
dir.create(plotPath, showWarnings = FALSE)

orga_type = tolower(organism)

#mouse code
if(orga_type == "mouse"){
  universe_mouse_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene_universe, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
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


#rn =4
cluster_meta_annot_excel =   function(rn, clusterdf) {
  dv= clusterdf[rn, ]
  print(paste0("Prediction for ", dv$color," containing ", dv$gene_no, " genes"))
  genes_dv = stri_split_regex(dv$gene_n, pattern=",") %>% unlist()

  ##mouse code
  if(orga_type == "mouse"){
    universe_orignal <- gene_universe
    genes_cluster = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unlist(genes_dv), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    cluster_genes <- genes_cluster[,2]
    list_of_genes_mouse_human <- list(cluster_genes)
    entrez_de = clusterProfiler::bitr(unlist(genes_dv), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    module_entrez_mouse <- unlist(entrez_de[2],use.names = FALSE)
    entrez_de = clusterProfiler::bitr(unlist(list_of_genes_mouse_human), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    module_entrez_mouse_human <- unlist(entrez_de[2],use.names = FALSE)


  }else{

  #humans
  #@thomas shobhit some symbols do not have entrezids, maybe they have aliases  e.g. C6orf48
  #select entrez ids for the genes
  module_entrez= universe_Entrez %>%
    #drop_na(ENTREZID) %>%
    filter(SYMBOL%in%genes_dv) %>%
    pull(ENTREZID)

}

  #AnnotationDbi::mapIds(org.Hs.eg.db,genes_dv, 'ENTREZID', 'SYMBOL')

  wb = openxlsx::createWorkbook()

  pdf(paste0(originalwd,cutoff_wd,"/clusterProfiler/ClusterProfiler_", dv$color, ".pdf", sep=""), onefile = F,
      width = 15, height = 15)
  font_size=8

  plot.new()
  grid::grid.newpage()
  grid::pushViewport(viewport(layout=grid.layout(nrow=3,ncol=2)))



#GMT
  #mouse code
  if(orga_type == "mouse"){
    egmt <- clusterProfiler::enricher(unlist(list_of_genes_mouse_human), TERM2GENE=c1_hallmark_genes, universe = universe_mouse_human, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
  }else{
  egmt <- clusterProfiler::enricher(genes_dv,
                                    TERM2GENE=gmtfile_hallmarks,
                                    universe = universe_Entrez$SYMBOL,
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "none",
                                    qvalueCutoff = 1.0)
  }

  if(!is.null(egmt)){
    hallmark_plot <- clusterProfiler::dotplot(egmt, font.size = font_size, title = "Hallmark enrichment", orderBy="GeneRatio")

    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(hallmark_plot, newpage = FALSE)
    popViewport()

    sheet <- openxlsx::addWorksheet(wb, sheetName = "Hallmark")
    openxlsx::writeDataTable(wb, sheet, egmt@result)
  }


#KEGG
  #mouse
  if(orga_type == "mouse"){
    kegg_enrich <- enrichKEGG(gene = module_entrez_mouse, organism = 'mmu', pvalueCutoff=0.05,universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)

    Enriched_Kegg<-NULL
    Enriched_Kegg_obj<-NULL
    reg_Mm=AnnotationDbi::select(org.Mm.eg.db,as.character(unlist(genes_dv)),"ENTREZID","SYMBOL",multiVals="first")

    if(!is.null(kegg_enrich))
    {
      if(nrow(summary(kegg_enrich))>0)
      {
        df_kk<-as.data.frame(summary(kegg_enrich))[1:8]
        for(x in 1:length(df_kk[,8]))
        {
          temp<-strsplit(df_kk[x,8],"/")
          id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
          df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
        }
        Enriched_Kegg<-df_kk
        Enriched_Kegg_obj<-kegg_enrich
      }
      else
      {
        Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 8))
        Enriched_Kegg_obj<-NULL
      }
    }

    if(!is.null(kegg_enrich) && kegg_enrich@result %>%nrow()>0){


      kegg_enrich@result$symbols=
        lapply(1:nrow(kegg_enrich@result), function(x) paste(bitr(strsplit(kegg_enrich@result[x,8],"/")[[1]], fromType="ENTREZID",toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL,collapse = "/"))

        # paste(bitr(strsplit(right@result[i,8],"/")[[1]], fromType="ENTREZID",toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL,collapse = "/")


      kegg_plot= clusterProfiler::dotplot(kegg_enrich,
                                          font.size= font_size,
                                          title="KEGG_enrichment",
                                          orderBy="GeneRatio")
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      print(kegg_plot, newpage= F)
      popViewport()
      sheet =  openxlsx::addWorksheet(wb, sheetName = "KEGG")
      openxlsx::writeDataTable(wb, sheet, kegg_enrich@result)

      #Enriched_Kegg = kegg_enrich@result
      #Enriched_Kegg_obj = kegg_enrich
    }

  }else{

  #human
  kegg_enrich = clusterProfiler::enrichKEGG(gene = module_entrez, organism = 'hsa',
                           pvalueCutoff=0.05,
                           universe = universe_Entrez$ENTREZID,
                           pAdjustMethod = "none",
                           qvalueCutoff = 1.0)

  Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 10))
  Enriched_Kegg_obj<-NULL

  if(!is.null(kegg_enrich) && kegg_enrich@result %>%nrow()>0){
    kegg_enrich@result$symbols=
      lapply(1:nrow(kegg_enrich@result),
             function(x) universe_Entrez %>%
               filter(ENTREZID%in%(stri_split_fixed(kegg_enrich@result[x,"geneID"], pattern = "/") %>% unlist())) %>%
               pull(SYMBOL) %>%paste0(collapse = "/"))  %>%
      unlist()
    kegg_plot= clusterProfiler::dotplot(kegg_enrich,
                                        font.size= font_size,
                                        title="KEGG_enrichment",
                                        orderBy="GeneRatio")
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(kegg_plot, newpage= F)
    popViewport()
    sheet =  openxlsx::addWorksheet(wb, sheetName = "KEGG")
    openxlsx::writeDataTable(wb, sheet, kegg_enrich@result)

    #Enriched_Kegg = kegg_enrich@result
    #Enriched_Kegg_obj = kegg_enrich
  }

}




#Reactome
  #mouse code
  if(orga_type == "mouse"){
    reactome_enrich <- ReactomePA::enrichPathway(gene=module_entrez_mouse, organism = "mouse", pvalueCutoff=0.05, universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)

    Enriched_ReacTome<-NULL
    Enriched_ReacTome_obj<-NULL

    if(!is.null(reactome_enrich))
    {
      if(nrow(summary(reactome_enrich))>0)
      {
        df_kk<-as.data.frame(summary(reactome_enrich))[1:8]
        for(x in 1:length(df_kk[,8]))
        {
          temp<-strsplit(df_kk[x,8],"/")
          id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
          df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
        }
        Enriched_ReacTome<-df_kk
        Enriched_ReacTome_obj<-reactome_enrich@result
      }
      else
      {
        Enriched_ReacTome<-data.frame(matrix(NA, nrow = 0, ncol = 8))
        Enriched_ReacTome_obj<-NULL
      }
      }

    reactome_enrich@result$symbols=
      lapply(1:nrow(reactome_enrich@result), function(x) paste(bitr(strsplit(reactome_enrich@result[x,8],"/")[[1]], fromType="ENTREZID",toType="SYMBOL", OrgDb="org.Mm.eg.db")$SYMBOL,collapse = "/"))


     }else{

  #human

  reactome_enrich <- ReactomePA::enrichPathway(gene=module_entrez, organism = "human",
                                        pvalueCutoff=0.05,
                                        universe = universe_Entrez$ENTREZID,
                                        #%>% drop_na(ENTREZID)%>%pull(ENTREZID),
                                        pAdjustMethod = "none",
                                        qvalueCutoff = 1.0)

  if(!is.null(reactome_enrich) & reactome_enrich@result %>%nrow()>0 )
  {
    reactome_enrich@result$symbols=
      lapply(1:nrow(reactome_enrich@result),
             function(x) universe_Entrez %>%
               filter(ENTREZID%in%(stri_split_fixed(reactome_enrich@result[x,"geneID"], pattern = "/") %>% unlist())) %>%
               pull(SYMBOL) %>%paste0(collapse = "/"))  %>%
      unlist()
  }


  }
  if(!is.null(reactome_enrich)){
    if(reactome_enrich@result %>%nrow()>0){
    # reactome_enrich@result$symbols=
    #   lapply(1:nrow(reactome_enrich@result),
    #          function(x) universe_Entrez %>%
    #            filter(ENTREZID%in%(stri_split_fixed(reactome_enrich@result[x,"geneID"], pattern = "/") %>% unlist())) %>%
    #            pull(SYMBOL) %>%paste0(collapse = "/"))  %>%
    #   unlist()
    reactome_plot= clusterProfiler::dotplot(reactome_enrich,
                                        font.size= font_size,
                                        title="reactome_enrichment",
                                        orderBy="GeneRatio")
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(reactome_plot, newpage= F)
    popViewport()
    sheet =  openxlsx::addWorksheet(wb, sheetName = "Reactome")
    openxlsx::writeDataTable(wb, sheet, reactome_enrich@result)

    }
  }


  ##DO


  if(orga_type == "mouse"){

    DO_enrich <- DOSE::enrichDO(gene=module_entrez_mouse_human,
                                ont = "DO",
                                pvalueCutoff=0.05,
                                universe = universe_Entrez_mouse_human,
                                pAdjustMethod = "none",
                                qvalueCutoff = 1.0)



  }else{


  DO_enrich <- DOSE::enrichDO(gene=module_entrez,
                              ont = "DO",
                              pvalueCutoff=0.05,
                              universe = universe_Entrez$ENTREZID,
                              pAdjustMethod = "none",
                              qvalueCutoff = 1.0)



  }

  if(!is.null(DO_enrich)){
    if(DO_enrich@result%>%nrow()>0){

    DO_enrich@result$symbols=
      lapply(1:nrow(DO_enrich@result), function(x) paste(bitr(strsplit(DO_enrich@result[x,8],"/")[[1]], fromType="ENTREZID",toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL,collapse = "/"))

    DO_plot= clusterProfiler::dotplot(DO_enrich,
                                      font.size= font_size,
                                      title="Disease enrichment",
                                      orderBy="GeneRatio")
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    print(DO_plot, newpage= F)
    popViewport()
    sheet =  openxlsx::addWorksheet(wb, sheetName = "Disease")
    openxlsx::writeDataTable(wb, sheet, DO_enrich@result)

  }}

#enrich GO


  if(orga_type == "mouse"){


    GO_enrich <- clusterProfiler::enrichGO(gene=module_entrez_mouse,
                                           ont = "BP",
                                           OrgDb = org.Mm.eg.db,
                                           pvalueCutoff=0.05,
                                           universe = universe_Entrez_mouse,
                                           pAdjustMethod = "none",
                                           qvalueCutoff = 1.0,
                                           readable = T)


  }else{

  GO_enrich <- clusterProfiler::enrichGO(gene=module_entrez,
                                         ont = "BP",
                                         OrgDb = org.Hs.eg.db,
                                         pvalueCutoff=0.05,
                                         universe = universe_Entrez$ENTREZID,
                                         pAdjustMethod = "none",
                                         qvalueCutoff = 1.0,
                                         readable = T)




  }
  if(!is.null(GO_enrich))
    #@shobhit add a column for gene symbol, shift the geneid in the go result to symbol and find the corresponding ENTREZID for geneId
  {
    enrichGO_plot <- clusterProfiler::dotplot(GO_enrich,
                                              font.size = font_size,
                                              title = "GO enrichment",
                                              showCategory=20,orderBy="GeneRatio"
                                              )
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(enrichGO_plot, new = FALSE)
    popViewport()

    sheet <- openxlsx::addWorksheet(wb, sheetName = "GO")
    openxlsx::writeDataTable(wb, sheet, GO_enrich@result)

  }


#TF prediction
  #weird
  #TF prediction
  # plotFunction<-function(){
  #   p<- ggplot2::ggplot(TFs_matrix_all_plot.melt.mean,aes(x=variable,y=mean,fill=merged)) +
  #     geom_bar(stat="identity",position = "dodge") +
  #     facet_wrap(~variable,scales = "free")
  #   print(p)
  # }

  if(orga_type == "mouse"){
    TFtable <- pcaGoPromoter::primo(genes_dv, inputType = "geneSymbol", org = "Mm")


  }else{


    TFtable <- pcaGoPromoter::primo(genes_dv, inputType = "geneSymbol", org = "Hs")


  }


  #TFtable$overRepresented$gene %>% as.character() %>% stri_replace_all_regex(pattern = "_.*", replacement = "")
  TFs=TFtable$overRepresented$gene %>%
    as.character() %>%
    stri_split_regex(pattern = "_") %>%
    #map(1) %>%
    unlist() %>%
    stri_split_regex(pattern = "::") %>%
    unlist() %>%
    unique()

  TFS_intersect = NULL


  capFirst <- function(s) {
    s <- tolower(s)
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
  }

 # macro_genes_Kevin$ont <- capFirst(macro_genes_Kevin$ont)

 if(!is.null(TFs) & length(TFs) >0){
   sheet <- openxlsx::addWorksheet(wb, sheetName = "TFoverrepresented")
   openxlsx::writeDataTable(wb, sheet, TFtable$overRepresented)
   TFS_intersect <- intersect(toupper(TFs),gene_universe)
    TFs_matrix_all_plot <- normdata[,c("merged",TFS_intersect)]
    TFs_matrix_all_plot.melt <- tidyr::gather(TFs_matrix_all_plot,gene, value, -merged )


    # calculate means
    # TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt %>%
    #   dplyr::group_by(merged,gene) %>%
    #   dplyr::summarise(mean=mean(value))

    ####plotFunction####
    plotTheThing<-function(){

      p<- ggplot(TFs_matrix_all_plot.melt,aes(x=merged,y=value,fill=merged)) +
        geom_boxplot(outlier.colour="red", outlier.shape=1, outlier.size=2)+
        geom_jitter(shape=16, position=position_jitter(0.2))+
        facet_wrap(~gene,scales = "free")

      # ggplot(TFs_matrix_all_plot.melt.mean,aes(x=gene,y=mean,fill=merged)) +
      #   geom_bar(stat="identity",position = "dodge") +
      #   labs(fill = paste0(group, collapse = "-"))+
      #   facet_wrap(~gene,scales = "free")



      print(p)

    }

    # here you put the plot into the excel sheet
    png(paste0(working_directory,cutoff_wd, "/clusterProfiler/",dv$color,".png"), width=1450, height = 450)
    plotTheThing()
    dev.off()
    openxlsx::insertImage(wb, sheet, paste0(working_directory,cutoff_wd, "/clusterProfiler/",dv$color,".png"))
  }


#  @thomas some weird pdf code

   if(length(TFS_intersect) != 0L){
     #@thomas why first 4 genes?
    TFs_matrix <- normdata[,c("merged",head(TFS_intersect,4))]
     #TFs_matrix.melt <- melt(TFs_matrix, id="merged")        # calculate means

     TFs_matrix.melt <- tidyr::gather_(TFs_matrix, "variable", "value",gather_col = c(head(TFS_intersect,4)) )
     TFs_matrix.melt %>% group_by(merged,variable)  -> TFs_matrix.melt.mean
     TFs_matrix.melt.mean$value<-as.numeric(sub(",", ".", TFs_matrix.melt.mean$value, fixed = TRUE))
     condis <- unique(normdata$merged)
     TFs_matrix.melt.mean$merged <- factor(TFs_matrix.melt.mean$merged,levels = condis)

  #   # for pdf
    TF_plot <- ggplot(TFs_matrix.melt.mean,aes(x=merged,y=value,fill=merged)) +
      geom_boxplot(outlier.colour="red", outlier.shape=1, outlier.size=2)+
      geom_jitter(shape=16, position=position_jitter(0.2))+
      facet_wrap(~variable,scales = "free")
  #
  #   # embed in pdf
     pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
     print(TF_plot, newpage = FALSE)
     popViewport()
   }
  dev.off()
  openxlsx::saveWorkbook(wb, paste0(originalwd,cutoff_wd,"/clusterProfiler/ClusterProfiler_",dv$color,".xlsx"), overwrite = T)

  }


lapply(1:nrow(cluster_data), cluster_meta_annot_excel, clusterdf=cluster_data)
#setwd(originalwd)

}
