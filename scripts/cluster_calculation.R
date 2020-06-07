cluster_calculation <- function(igraph,
                                cluster_algo,
                                no_of_iterations,
                                max_cluster_count_per_gene,
                                min_cluster_size,
                                GFC) {


  g = igraph

  comps = count_components(g)

  cluster_algo_list =c("cluster_label_prop",
                       "cluster_fast_greedy",
                       "cluster_louvain",
                       "cluster_infomap",
                       "cluster_walktrap")

  algos_to_use = switch(cluster_algo=="auto", cluster_algo_list, cluster_algo)


  ##choosing algorithm for clustering
  ##function to calculate modularity score
  #@thomas: point of iteration and seeding

  seed=.Random.seed

  # best_algo_df = data.frame(best_algo="", mod_score=0)
  # gwc = NULL

  cluster_calculations =function(graph_obj, algo, case, iter) {
    ###using function name (string) as a function
    seeda=seed[iter]
    set.seed(seeda)
    print(seeda)
    cfg= get(algo)(graph_obj)
    #print(cfg$names)
    mod_score =modularity(graph_obj, cfg$membership)
    mod_df= data.frame(modularity_score=mod_score, cluster_algorithm=algo, stringsAsFactors = F)

    ##making switch so that in the end when only the best algorithm is to be used then the same function can be used
    output = switch(case, best= cfg$membership, test= mod_df, final=cfg)

    print(paste0(algo,"algorithm tested"))
    return(output)
  }

  if(cluster_algo=="auto"){

  df_modularity_score = do.call("rbind", lapply(algos_to_use,
                                                cluster_calculations,
                                                graph_obj=g,
                                                case="test", iter=1))


  cluster_algo_used= df_modularity_score %>%
    dplyr::filter(modularity_score==max(modularity_score)) %>%
    dplyr::select(cluster_algorithm) %>%
    as.character()

  print(paste(cluster_algo_used, "will be used based on your input (if not auto option was specified) or the highest modularity score "))

  }else{
    cluster_algo_used <- cluster_algo
    print(paste(cluster_algo_used, "will be used based on your input (if not auto option was specified) or the highest modularity score "))

  }


  igraph_list = list()
  igraph_list[[1]] = g

  ###apply the best clustering algorithm
  gene_which_cluster=do.call("cbind", lapply(1:no_of_iterations,
                                             cluster_calculations,
                                             algo=cluster_algo_used,
                                             case="best",
                                             graph_obj=g))



  ##frequency and identity of cluster assingment of genes
  if(base::ncol(gene_which_cluster) > 1) {
    gene_cluster_ident = apply(gene_which_cluster,1, function(x){
      if(length(unique(x)) > max_cluster_count_per_gene) {   #LISA: was >=
        0
      }else{
        names(which(table(x) == max(table(x))))[1]
      }
    })
  } else{ gene_cluster_ident = gene_which_cluster[,1]}
  #gene_which_cluster <- as.matrix(gene_which_cluster)

  white_genes_clustercounts <- as.integer(grep(gene_cluster_ident, pattern = "\\b0\\b") %>%  #LISA: added \\b to avoid that also 20,... are grepped
    length() %>% as.character())


  print(paste(white_genes_clustercounts, "genes were assigned to more than", max_cluster_count_per_gene,
              "clusters. These genes are assigned to Cluster 0 and will be painted white in the network."))



  cluster_Data = data.frame(genes=vertex_attr(g, "name"),
                            clusters= paste0("Cluster ",gene_cluster_ident),
                            stringsAsFactors = FALSE)

  #summarize the data
  #produces a table where col are cluster name, number of components,
  #names of genes in cluster

  dfk=cluster_Data %>%
    dplyr::count(clusters,genes) %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(gene_no= sum(n), gene_n = paste0(genes,collapse = ",")) %>%
    dplyr::mutate(cluster_included=ifelse(gene_no>=min_cluster_size, "yes", "no"), color="white")



  #LISA: different color options (replacement for WGCNA)

  ##ggplot

  color.cluster <- c("orchid", "maroon", "darkgreen",  "darkorange", "darkgrey", "gold", "steelblue", "indianred",
                     "pink", "lightgreen", "lightblue","sandybrown",   "khaki",  "turquoise","darkblue",
                     "cadetblue","greenyellow","cyan", "thistle", "darkmagenta", "coral", "red", "blue",
                     "green", "yellow", "brown", "black")

  plot_clusters <- ggplot(data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ],
                          aes(x = clusters)) +
    geom_bar(aes(fill = clusters)) +
    scale_fill_manual(values = color.cluster)

  plot_clust <- ggplot_build(plot_clusters)

  dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color" ] <- plot_clust$data[[1]]["fill"]


  ### cols25 palette and modulecolor package
  # dfk$clus_number <- dfk$clusters
  # dfk$clus_number <- as.numeric(gsub("[a-zA-Z]", "", dfk$clus_number))
  # dfk[dfk$included=="yes","color"] =
  #   moduleColor::labels2colors(dfk$clus_number[dfk$included=="yes"],
  #                              colorSeq = pals::cols25())
  # dfk$clus_number <- NULL


  ###standard colors and modulecolor package
  # dfk$clus_number <- dfk$clusters
  # dfk$clus_number <- as.numeric(gsub("[a-zA-Z]", "", dfk$clus_number))
  # dfk[dfk$included=="yes","color"] =
  #   moduleColor::labels2colors(dfk$clus_number[dfk$included=="yes"])




  white_genes_clustersize <- as.integer(dfk %>% dplyr::filter(cluster_included=="no")%>%
                                          dplyr::summarise(n=sum(gene_no)) %>% purrr::map(1))


  print(paste0(white_genes_clustersize, " genes were assigned to clusters with a smaller size than the defined minimal cluster size of ",
             min_cluster_size, " genes per cluster. These genes will also be painted white in the network."))




  # for each cluster produces a row of means per condition (info data) for all genes within the cluster
  #included_clusters = subset(dfk, included=="yes")
  cluster_df=dfk
  gfc_dat = GFC

  gfc_mean_clustergene=function(rownum, cluster_df, gfc_dat){
    d1 = cluster_df[rownum,]
    gene_names= d1["gene_n"] %>%
      stri_split_regex(pattern = ",") %>%
      unlist()
    gfc_means = gfc_dat[gfc_dat$Gene%in%gene_names,] %>%
      dplyr::select(-Gene) %>%
      colMeans()
    d1$conditions = paste0(names(gfc_means), collapse = "#")
    d1$grp_means = paste0(round(gfc_means,3) , collapse = ",")
    return(d1)
  }


  dfk_allinfo=do.call("rbind", lapply(1:nrow(dfk), gfc_mean_clustergene, cluster_df= dfk,
                                      gfc_dat=GFC))
  dfk_allinfo$vertexsize = ifelse(dfk_allinfo$cluster_included=="yes",3,1)



  return(dfk_allinfo)

}

cluster_info_long_format <- function(cluster_information) {

    s <- strsplit(cluster_information$gene_n, split = ",")
    cluster_info_long_format <- data.frame(color = rep(cluster_information$color, sapply(s, length)), gene_n = unlist(s))
    return(cluster_info_long_format)

}

