##collect function parameters
allargs <- function(orig_values = FALSE) {
  # get formals for parent function
  parent_formals <- formals(sys.function(sys.parent(n = 1)))
  
  # Get names of implied arguments
  fnames <- names(parent_formals)
  
  # Get currently set values for named variables in the parent frame
  args <- evalq(as.list(environment()), envir = parent.frame())
  return(args)
}
#############To flatten Correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


###cor funktioniert vllt erst auf pvalues cutten oder 0.5 correlation

#############Heatmap
# Heatmap_Lea<-function(x,filter){
#   if(is.numeric(x)==T){
#     rownames(Top_genes)<-Top_genes$ID
#     Top_genes$ID<-NULL
#     Heatmap_data<-t(Top_genes[,(x-1):ncol(Top_genes)])
#   }else{
#     print(c("x must be a number, corresponds to column where first expression data is in"))
#   }
#   if(is.character(filter)==T){
#     Heatmap_data_filtered<-Heatmap_data[filter %in% rownames(Heatmap_data),]
#     aheatmap(Heatmap_data_filtered,scale="row", color = "-RdBu", annColors = "Set3", 
#              annCol =Top_genes$merged, 
#              Rowv = NA )
#   }else{
#     print(c("filter must be characters (e.g. list of genes)"))
#   }
# }

##r²-values for scale free topology
r_squared_values = function(graph) {# calculate degree
  d = igraph::degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  if(length(probability)==0){
    R.square<-0
  }else{
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
  }
  # cozf = coef(reg)
  # power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  # alpha = -cozf[[2]]
  # if(summary(reg)$r.squared==c("NaN")){
  #   R.square<-0
  # }else{
  # R.square = summary(reg)$r.squared
  # }
  
  #print(paste("Alpha =", round(alpha, 3)))
  #print(paste("R square =", round(R.square, 3)))
  r_square_val<-R.square}

fit_power_law = function(graph,pdf=T) {
  # calculate degree
  d = igraph::degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  if(pdf){
    pdf("Degree_distribution_plot.pdf")
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
  dev.off()
  }else{
    plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
         col = 1, main = "Degree Distribution")
    curve(power.law.fit, col = "red", add = T, n = length(d))
  }
}

###change vertex.frame.width by creating new shape
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

#igraph needs to be loaded beforehand
add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                                vertex.frame.width=1))

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}

scale_free_top = function(graph) {
  # calculate degree
  d = igraph::degree(graph, mode = "all")
  dd = igraph::degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  R.square = round(R.square, 3)
  
  data.frame(R.squared=R.square,degree=degree,Probs=probability)
  #print(paste("Alpha =", round(alpha, 3)))
  #print(paste("R square =", round(R.square, 3)))
  # plot
  #plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
  #col = 1, main = "Degree Distribution")
  #curve(power.law.fit, col = "red", add = T, n = length(d))
}


######correlation
correlation<-function(type_of_correlation=c("pearson","spearman"),pVal_treshold=0.05,save_raw_data=FALSE){
  #summary[[c("correlation")]]<-data.frame(type_of_correlation,pVal_treshold,save_raw_data)
  summary<-allargs()
  original_data_t <- t(original_data)
  cor_mat_norm_data <- as.matrix(original_data_t)
  ##funktioniert mit großem datensatz! flattencormatrix das problem?
  if(type_of_correlation==c("pearson")){
    cor_val_norm_data <- Hmisc::rcorr(cor_mat_norm_data , type = "pearson")
  }
  
  if(type_of_correlation==c("spearman")){
    cor_val_norm_data <- rcorr(cor_mat_norm_data , type = "spearman")
  }
  
  pvalues_cor_norm_data <- cor_val_norm_data$P
  if(is.function(flattenCorrMatrix)==TRUE){
    flat_cor_mat <- flattenCorrMatrix(cor_val_norm_data$r , cor_val_norm_data$P)
    flat_cor_mat_sig <- flat_cor_mat[flat_cor_mat$p < pVal_treshold , ]
    modul_data_cor <- data.frame(from = flat_cor_mat_sig$row , to = flat_cor_mat_sig$column , correlation = flat_cor_mat_sig$cor , pval = flat_cor_mat_sig$p)
    rm(flat_cor_mat)
    rm(flat_cor_mat_sig)
    raw_data <- modul_data_cor
    raw_data <- raw_data[raw_data$correlation > 0 , ]
    #raw_data$group <- c("weight")
    if(save_raw_data==TRUE){
      title<-c("raw_data.txt")
      setwd(originalwd)
      write.table(raw_data , title , sep="\t")
      print(paste0("saved in ",getwd(), ". Filename:",title))
    }
  }else{
    print(" You must load the function flattenCorrMatrix!")
  }
  
  
  RT<-list()
  
  RT[["summary"]]<-summary
  RT[["raw_data"]]<-raw_data
  return(RT)
  
}

###corr
cutoff_visualisation<-function(correlation_df=correlation_df,min_corr=0.5,range_cutoff_length=10,min_no_for_cluster=10){
  #summary[[c("cutoff_visualisation")]]<-data.frame(min_corr,range_cutoff_length,min_no_for_cluster)
  summary<-allargs()
  
  correlation_df<-correlation_df[!(is.na(correlation_df$from)),]
  range_cutoff<-seq(from = min_corr , to = max(correlation_df$correlation) , length.out = range_cutoff_length)
  range_cutoff<-round(range_cutoff, 3)
  cut_off_data_frame<-data.frame(cut_off=range_cutoff,
                                 number_edges=1:length(range_cutoff),
                                 number_nodes=1:length(range_cutoff),
                                 r_squared_all=1:length(range_cutoff),
                                 #r_squared_biggest=1:length(range_cutoff),
                                 no_network=1:length(range_cutoff))
  
  
  data_for_scale_free_plot<-list()
  no<-0
  for(cutoff in range_cutoff){
    no<-no+1
    list_edges_nodes<-correlation_df[(correlation_df$correlation > cutoff) , ]
    if(nrow(list_edges_nodes)==0){
      cut_off_data_frame$r_squared_all[no]<-0
      cut_off_data_frame$no_network[no]<-0
      cut_off_data_frame$number_edges[no]<-0
      cut_off_data_frame$number_nodes[no]<-0
      data_for_scale_free_plot[[no]]<-data.frame(R.squared=0,
                                                 degree=0,
                                                 Probs=0,
                                                 cutoff=cutoff,
                                                 no_edges=cut_off_data_frame$number_edges[no],
                                                 no_nodes=cut_off_data_frame$number_nodes[no],
                                                 no_of_networks=cut_off_data_frame$no_network[no])
    }else{
      igraph_object<-graph_from_data_frame(list_edges_nodes,vertices = NULL,directed=FALSE)
      cut_off_data_frame$r_squared_all[no] <- r_squared_values(igraph_object)
      graph_components <- clusters(igraph_object)
      cut_off_data_frame$no_network[no]<-length(graph_components$csize[graph_components$csize>min_no_for_cluster])
      cut_off_data_frame$number_edges[no]<-length(list_edges_nodes$from)
      cut_off_data_frame$number_nodes[no]<-length(unique(c(as.character(list_edges_nodes[,1]),as.character(list_edges_nodes[,2]))))
      data_for_scale_free_plot[[no]]<-data.frame(scale_free_top(igraph_object),
                                                 cutoff=cutoff,
                                                 no_edges=cut_off_data_frame$number_edges[no],
                                                 no_nodes=cut_off_data_frame$number_nodes[no],
                                                 no_of_networks=cut_off_data_frame$no_network[no])
    }
    
    print(paste0(no," out off ",range_cutoff_length, "cutoffs are calculated" ))
  }
  data_for_scale_free_plot<-list.rbind(data_for_scale_free_plot)
  library(plotly)
  library(crosstalk)
  g<- crosstalk::SharedData$new(data_for_scale_free_plot,~cutoff)
  p1.0<-plot_ly(data_for_scale_free_plot, 
                y=~no_nodes,
                x=~cutoff,
                #type="scatter",
                #mode="scatter",
                hoverinfo= "y")%>%
    add_markers(alpha=0.1, color=I("black"))%>%
    add_markers(data=g, frame=~cutoff,color=I("red"))%>%
    plotly::layout(yaxis = list(zeroline=FALSE,title="No. of Nodes"),
                   xaxis = list(zeroline = FALSE),range= c(min_corr,1))


  #,xaxis = list(zeroline = FALSE,range= c(min_corr,1)))
  
  p1.1<-plot_ly(data_for_scale_free_plot, 
                y=~no_edges,
                x=~cutoff,
                #type="scatter",
                #mode="scatter",
                hoverinfo = "y")%>%
    add_markers(alpha=0.1, color=I("black"))%>%
    add_markers(data=g, frame=~cutoff, color=I("red"))%>%
    plotly::layout(yaxis = list(zeroline = FALSE,title=c("No. of Edges")),
           xaxis = list(zeroline = FALSE),
           range= c(min_corr,1))
  p1.2<-plot_ly(data_for_scale_free_plot,
                y=~no_of_networks,
                x=~cutoff,
                #type="scatter",
                #mode="scatter",
                hoverinfo = "y") %>%
    add_markers(alpha=0.1, color=I("black"))%>%
    add_markers(data=g, frame=~cutoff, color=I("red"))%>%
    plotly::layout(yaxis = list(zeroline = FALSE,
                        title=c("No. of Networks")),
           xaxis = list(zeroline = FALSE),
           range= c(min_corr,1))
  
  p1.3<-plot_ly(data_for_scale_free_plot,
                y=~R.squared,
                x=~cutoff,
                #type="scatter",
                #mode="scatter",
                hoverinfo = "y") %>%
    add_markers(alpha=0.1, color=I("black"))%>%
    add_markers(data=g, frame=~cutoff, color=I("red"))%>%
    plotly::layout(yaxis = list(zeroline = FALSE,
                        title="r² Value"),
           xaxis = list(zeroline = FALSE,
                        range= c(min_corr,1)))
  
  p2<-plot_ly(data_for_scale_free_plot,
              x=~degree,
              y=~Probs,
              text=data_for_scale_free_plot$R.squared,
              hoverinfo ="text") %>%
    add_markers(color=I("black"),alpha=0.1) %>%
    add_markers(data=g, frame=~cutoff, color=I("red"))%>%
    plotly::layout(xaxis=list(type="log"),yaxis=list(type="log"))
  
  p_1<-plotly::subplot(p1.3,p2, widths  = c(0.3, 0.7), titleY = TRUE,titleX = TRUE) 
  p_2<-plotly::subplot(p1.2,p1.1,p1.0, nrows = 1, widths = c(0.2, 0.4,0.4),
               titleY = TRUE,titleX = TRUE) 
  
  end_plot<-plotly::subplot(p_1,p_2, nrows = 2,heights = c(0.7,0.3), titleY = TRUE,
                    titleX = TRUE) %>%
    hide_legend() %>%
    animation_opts(frame=100,transition=0,easing="linear", redraw = FALSE) %>%
    #layout(hovermode = "y", margin = list(l = 100)) %>%
    #highlight("plotly_selected", color = "blue", opacityDim = 1, hoverinfo = "none")%>%
    animation_slider(
      currentvalue = list(prefix = "Cutoff ", font = list(color="red"))
    )
  list_cutoff_options<-list()
  list_cutoff_options[[c("Plotly_object")]]<-end_plot
  list_cutoff_options[[c("Cutoff_df")]]<-cut_off_data_frame
  list_cutoff_options[["summary"]]<-summary
  return(list_cutoff_options) 
}

#optimal cutoff
optimal_cutoff<-function(cutoff_table=cutoff_options){
  
  
  n <- length(unique(cutoff_table$Cutoff_df$no_network))
  num_of_subnetworks <- sort(unique(cutoff_table$Cutoff_df$no_network))[2]
  
  optimal_cutoff <- cutoff_table$Cutoff_df %>% 
    filter(no_network == num_of_subnetworks)

  optimal_cutoff <- optimal_cutoff %>%
    filter(number_edges < 30000)
    
    optimal_cutoff <- optimal_cutoff %>%
    filter(number_nodes > 0.3*max(optimal_cutoff$number_nodes)) %>%
    arrange(desc(r_squared_all)) %>%
    filter(row_number()==1)
  
  return(optimal_cutoff$cut_off) 
}

##plotnetwork
plot_network<-function(data=correlation_df,layout=layout_with_kk, min_nodes_number_for_network=10,show_HC= TRUE, pdf = T){
  summary<-allargs()
  #how include layout in summary?
  #summary[[c("plot_network")]]<-data.frame(min_nodes_number_for_network,show_HC)
  cutted_data<-data[data$correlation>chosen_cutoff,]
  igraph_object<-igraph::graph_from_data_frame(cutted_data,directed=FALSE)
  
  if(show_HC==TRUE){
    if(pdf){
      
      nodes_list<-V(igraph_object)
      tmp<-original_data[nodes_list,]
      tmp <- tmp[apply(tmp, MARGIN = 1, FUN = function(x) sd(x) != 0),]
      tmp <- tmp[complete.cases(tmp),]
      
      pdf("HC_network_genes.pdf")
      
      pheatmap::pheatmap(tmp, 
                         cluster_row = T,
                         cluster_cols = T,
                         color = col.pal,
                         scale = c("row"),
                         annotation_col = info_Dataset)
      dev.off()
    }else{
      nodes_list<-V(igraph_object)
      tmp<-original_data[nodes_list,]
      tmp <- tmp[apply(tmp, MARGIN = 1, FUN = function(x) sd(x) != 0),]
      tmp <- tmp[complete.cases(tmp),]
      
      pdf("HC_network_genes.pdf")
      
      pheatmap::pheatmap(tmp, 
                         cluster_row = T,
                         cluster_cols = T,
                         color = col.pal,
                         scale = c("row"),
                         annotation_col = info_Dataset)
      dev.off()
      
    }
  }
  if(pdf){
    if(count_components(igraph_object)==1){
      pdf("Selected_network_layout.pdf")
      igraph_list<-list()
      igraph_list[[1]]<-igraph_object
      layout_df<-BiocGenerics::do.call(layout, igraph_list)
    
      g<-igraph_object
      plot(igraph_object,layout=layout_df,
           vertex.size = 3 ,
          vertex.frame.width = 3/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey")
      dev.off()
      
    }else{
      components<-1:igraph::components(igraph_object)[[3]]
      components_numbers<-components[igraph::components(igraph_object)[[2]]>min_nodes_number_for_network]
      V(igraph_object)$comp<-igraph::components(igraph_object)[[1]]
      list_of_components_graphs<-list()
      for(component in components_numbers){
        list_of_components_graphs[[as.character(component)]] <- induced_subgraph(igraph_object , V(igraph_object)$comp == component)
      }
      layouts_on_list <- lapply(list_of_components_graphs , layout)
      layout_df <- merge_coords(list_of_components_graphs , layouts_on_list)
      g <- disjoint_union(list_of_components_graphs)
      
      pdf("Selected_network_layout.pdf")
      plot(g,layout=layout_df,
           vertex.size = 3 ,
           vertex.frame.width = 3/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey") 
      dev.off()
      
    }
  }else{
    if(count_components(igraph_object)==1){
      pdf("Selected_network_layout.pdf")
      igraph_list<-list()
      igraph_list[[1]]<-igraph_object
      layout_df<-BiocGenerics::do.call(layout, igraph_list)
      
      g<-igraph_object
      plot(igraph_object,layout=layout_df,
           vertex.size = 3 ,
           vertex.frame.width = 3/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey")

    }else{
      components<-1:igraph::components(igraph_object)[[3]]
      components_numbers<-components[igraph::components(igraph_object)[[2]]>min_nodes_number_for_network]
      V(igraph_object)$comp<-igraph::components(igraph_object)[[1]]
      list_of_components_graphs<-list()
      for(component in components_numbers){
        list_of_components_graphs[[as.character(component)]] <- induced_subgraph(igraph_object , V(igraph_object)$comp == component)
      }
      layouts_on_list <- lapply(list_of_components_graphs , layout)
      layout_df <- merge_coords(list_of_components_graphs , layouts_on_list)
      g <- disjoint_union(list_of_components_graphs)
      
      plot(g,layout=layout_df,
           vertex.size = 3 ,
           vertex.frame.width = 3/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey") 

    }
  }
  
  return_list<-list()
  return_list[[c("graph_object")]]<-g
  return_list[[c("layout")]]<-layout_df
  return_list[[c("summary")]]<-summary
  return(return_list)
}

#cluster_algo<-c("cluster_label_prop")
heatmap_clustered<-function(igraph_object=igraph_object,
                            cluster_algo=c("auto","cluster_label_prop"),
                            layout_for_network=layout,
                            iterations=TRUE,
                            no_of_iterations=10,
                            max_cluster_count_per_gene=8,
                            min_cluster_size=10,
                            desicion_to_see_plot=TRUE,
                            print_to_pdf=FALSE,
                            name_for_pdf=c("Heatmap_greater_clusters_clustered"),
                            average = "mean"){
  # summary[[c("heatmap_clustered")]]<-data.frame(iterations,no_of_iterations,max_cluster_count_per_gene,
  #                                      min_cluster_size,desicion_to_see_plot,desicion_to_save_plot,name_for_pdf_plot,
  #                                      print_to_pdf,name_for_pdf)
  summary<-allargs()
  g<-igraph_object
  if(cluster_algo==c("auto") & components(g)$no==1){
    df_modularity_score<-data.frame(modularity_score=1:5,cluster_algo=c("cluster_label_prop","cluster_fast_greedy","cluster_louvain","cluster_infomap",
                                                                        "cluster_walktrap"))
    cfg <- cluster_label_prop(g)
    df_modularity_score$modularity_score[1] <-modularity(g,cfg$membership)
    print("1/5 cluster algorithms tested")
    cfg <- cluster_fast_greedy(g)
    df_modularity_score$modularity_score[2] <-modularity(g,cfg$membership)
    print("2/5 cluster algorithms tested")
    cfg <- cluster_louvain(g)
    df_modularity_score$modularity_score[3] <-modularity(g,cfg$membership)
    print("3/5 cluster algorithms tested")
    cfg<-cluster_infomap(g)
    df_modularity_score$modularity_score[4] <-modularity(g,cfg$membership)
    print("4/5 cluster algorithms tested")
    cfg<-cluster_walktrap(g)
    df_modularity_score$modularity_score[5] <-modularity(g,cfg$membership)
    print("5/5 cluster algorithms tested")
    df_modularity_score<-df_modularity_score[order(df_modularity_score$modularity_score,decreasing = TRUE),]
    cluster_algo <-as.character(df_modularity_score$cluster_algo[1])
    print(paste0(cluster_algo, " was used to cluster, having the highest modularity score:",round(df_modularity_score$modularity_score[df_modularity_score$cluster_algo==cluster_algo],digits = 3)))
  }else{
    if(cluster_algo==c("auto") & components(g)$no>1){
      df_modularity_score<-data.frame(modularity_score=1:5,cluster_algo=c("cluster_label_prop","cluster_fast_greedy","cluster_louvain","cluster_infomap",
                                                                          "cluster_walktrap"))
      cfg <- cluster_label_prop(g)
      df_modularity_score$modularity_score[1] <-modularity(g,cfg$membership)
      print("1/5 cluster algorithms tested")
      cfg <- cluster_fast_greedy(g)
      df_modularity_score$modularity_score[2] <-modularity(g,cfg$membership)
      print("2/5 cluster algorithms tested")
      cfg <- cluster_louvain(g)
      df_modularity_score$modularity_score[3] <-modularity(g,cfg$membership)
      print("3/5 cluster algorithms tested")
      cfg<-cluster_infomap(g)
      df_modularity_score$modularity_score[4] <-modularity(g,cfg$membership)
      print("4/5 cluster algorithms tested")
      cfg<-cluster_walktrap(g)
      df_modularity_score$modularity_score[5] <-modularity(g,cfg$membership)
      print("5/5 cluster algorithms tested")
      
      df_modularity_score<-df_modularity_score[order(df_modularity_score$modularity_score,decreasing = TRUE),]
      cluster_algo <-as.character(df_modularity_score$cluster_algo[1])
      print(paste0(cluster_algo, " was used to cluster, having the highest modularity score:",round(df_modularity_score$modularity_score[df_modularity_score$cluster_algo==cluster_algo],digits = 3)))
      
    }
  }
  igraph_list<-list()
  igraph_list[[1]]<-g
  
  if(iterations==TRUE){
    gene_which_cluster<-list()
    for(rep in 1:no_of_iterations){
      if(rep==1){
        tic("time")
        
        clp<-BiocGenerics::do.call(cluster_algo, igraph_list)
        gene_which_cluster[[rep]]<-as.data.frame(clp$membership)
        time<-toc()
        print(paste0("Duration:",(time$toc-time$tic)*no_of_iterations,"secs approximatly"))
      }else{
        clp<-BiocGenerics::do.call(cluster_algo, igraph_list)
        gene_which_cluster[[rep]]<-as.data.frame(clp$membership)
        print(paste0(rep, "done out of",no_of_iterations))
      }
    }
    gene_which_cluster <- list.cbind(gene_which_cluster)
    
  }else{
    clp<-BiocGenerics::do.call(cluster_algo, igraph_list)
    gene_which_cluster[[1]]<-as.data.frame(clp$membership)
    gene_which_cluster <- list.cbind(gene_which_cluster)
  }
  
  if(iterations==TRUE){
    #gene_which_cluster <- as.matrix(gene_which_cluster)
    
    gene_which_cluster$cluster_max <- apply(gene_which_cluster , 1 , function(x){
      if(length(unique(x)) >= max_cluster_count_per_gene) {
        0
      }else{
        names(which(table(x) == max(table(x))))[1]
      }
    })
  }else{
    gene_which_cluster<-as.data.frame(gene_which_cluster[[1]])
    colnames(gene_which_cluster)<-c("cluster_max")
  }
  gene_which_cluster$cluster_max <- as.numeric(gene_which_cluster$cluster_max)
  for_printing<-as.data.frame(gene_which_cluster$cluster_max,check.names=FALSE)
  p<-as.data.frame(table(for_printing[,1] == 0))
  print(paste0(p[grepl("TRUE",p$Var1),2]," Genes are advised to more than ",max_cluster_count_per_gene," Clusters, and will be painted white in network-plot, summaraized in cluster 'waste' in heatmap and left out in upcoming analysis "))
  clp$membership <- gene_which_cluster$cluster_max  
  #gene_which_cluster$cluster_max <- as.numeric(gene_which_cluster$cluster_max)
  
  
  big_cluster <- list()
  number_left_out_due_tosmall_cluster<-list()
  cluster_min_sized_names <- list()
  big_cluster_names <- paste0("Cluster" ,unique(clp$membership))
  cluster_name_no<-0
  for(big_clust in unique(clp$membership)){
    cluster_name_no<-cluster_name_no+1
    big_clust<-as.numeric(big_clust)
    if (length(clp$membership[clp$membership==big_clust]) >= min_cluster_size) {                                            # here adjust the minimal size of a cluster that you still would like to include
      big_cluster[[big_cluster_names[cluster_name_no]]] <- as.data.frame(clp$names[clp$membership==big_clust])
      cluster_min_sized_names[[cluster_name_no]] <- big_clust
    }else{
      number_left_out_due_tosmall_cluster[[cluster_name_no]]<-length(clp$membership[clp$membership==big_clust])
    }
  }
  if(length(number_left_out_due_tosmall_cluster)>1){
    for_printing_2<-list.rbind(number_left_out_due_tosmall_cluster)
    sumx<-apply(for_printing_2,2,sum)
    print(paste0(sumx," Genes in total are painted white in network and left out in network due do their clustersize is < ",min_cluster_size))
  }
  big_cluster_all <- as.data.frame(list.rbind(big_cluster))
  cluster_min_sized_names <- as.data.frame(list.rbind(cluster_min_sized_names))
  
  list_of_means <- list()
  list_of_GFC_per_cluster<-list()
  names_of_GC_per_cluster <- paste0("Cluster" , unique(cluster_min_sized_names$V1))
  cluster_name_no<-0
  for(cluster in unique(cluster_min_sized_names$V1)) {
    cluster_name_no<-cluster_name_no+1
    #list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] <- GFC_all_genes[grepl(paste0(paste0("^" , clp[[cluster]] , "$") , collapse = "|") , GFC_all_genes$Gene) , ]
    list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]]<-GFC_all_genes[GFC_all_genes$Gene %in% (clp$names[clp$membership==cluster]),]
    list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]][ , ncol(GFC_all_genes)] <- NULL
    list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]] <- rbind(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]] , apply(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]] , 2 , average))
    list_of_means[[names_of_GC_per_cluster[cluster_name_no]]] <- apply(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster_name_no]]] , 2 , average)
  }
  list_of_means_new <- as.data.frame(list.rbind(list_of_means))
  
  #list_of_means_new <- list_of_means[cluster_min_sized_names$V1 , ]
  row.order <- hclust(dist(list_of_means_new))$order
  col.order <- hclust(dist(t(list_of_means_new)))$order
  dat_new <- list_of_means_new[row.order , col.order]
  
  df_molten <- reshape2::melt(as.matrix(dat_new))
  names(df_molten) <- c("Cluster" , "Tissue" , "GFC")
  df_molten$Cluster <- gsub(pattern="Cluster",replacement="",df_molten$Cluster)
  
  df_molten$Gen_No_cluster <- apply(df_molten , 1 , function(x){
    x <- as.numeric(x[1])
    length(clp$membership[clp$membership==x])
  })
  df_molten$Gen_No_cluster <- as.numeric(df_molten$Gen_No_cluster)
  df_molten$NEW_1 <- paste0(df_molten$Cluster , "[" , df_molten$Gen_No_cluster , "]")
  df_molten$Colors <- WGCNA::labels2colors(as.numeric(df_molten$Cluster) , zeroIsGrey = TRUE)
  
  if(c("grey") %in% unique(df_molten$Colors)){
    df_molten$Colors[grepl("^grey$",df_molten$Colors)]<-c("waste")
  }
  
  
  df_molten$NEW_2 <- paste0(df_molten$Colors , "[" , df_molten$Gen_No_cluster , "]")
  df_molten$Vertex_size <- as.numeric(c("3"))
  
  dat_new_dendo<-dat_new
  rownames(dat_new_dendo)<-unique(df_molten$NEW_2)
  
  
  
  #col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  # new_col = colorRamp2(breaks = c(min(dat_new_dendo), 0, max(dat_new_dendo)), colors = c("#053061", "#F7F7F7", "#67001F"))
  # if(print_to_pdf==TRUE){
  #   mainDir <- originalwd
  #   subDir <- chosen_cutoff
  #   dir.create(file.path(mainDir , subDir))
  #   setwd(file.path(mainDir , subDir))
  #   pdf(paste0(name_for_pdf,".pdf") , width = 8 , height = 8)
  #   draw(Heatmap(dat_new_dendo,
  #                #kmeans_k = 7,
  #                cluster_rows = T,
  #                cluster_columns = T,
  #                color = new_col))
  #   dev.off()
  #   print(paste0("PDF sent to ",getwd()))
  # }
  # 

  if(print_to_pdf==TRUE && desicion_to_see_plot ==TRUE){
    
    color_cluster_min_size <- as.data.frame(clp$membership)
    cluster_min_sized_names$V1 <- as.numeric(cluster_min_sized_names$V1)
    color_cluster_min_size$try <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        tmp <- df_molten[grep(paste0("^" , x[1] , "$") , df_molten$Cluster) , ]
        tmp[1 , 6]
      }else{
        "waste"
      }
    })
    cluster_min_sized_names$V1 <- as.character(cluster_min_sized_names$V1)
    color_cluster_min_size$`clp$membership` <- as.character(color_cluster_min_size$`clp$membership`)
    color_cluster_min_size$vertex_size <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        3
      }else{
        1
      }
    })
    hc_add_on <- color_cluster_min_size
    par(mfrow=c(1,2))
    frame()
    heatmap<-Heatmap(dat_new_dendo, 
                     cluster_rows = T,
                     cluster_columns = T,
                     color = new_col)
    gb = grid.grabExpr(draw(heatmap))#color auch hier ersetzt      pdf(name_for_netowrk_pdf , width = 8 , height = 8)
    
    pdf(paste0(name_for_pdf,".pdf") , width = 8 , height = 8)
    
    grid.arrange(gb,ncol=2)
    
    Sys.sleep(5)
    V(g)$size<-color_cluster_min_size$vertex_size
    color_cluster_min_size$vertex_color <- color_cluster_min_size$try
    color_cluster_min_size$vertex_color <- apply(color_cluster_min_size,1, function(x){
      
      if(x[4]=="waste"){"white"}else{x[4]}
      
    })
    #color_cluster_min_size[color_cluster_min_size$try%in%c("waste"),]<-c("grey")
    plot(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    
    #ppp<-ggnet2(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    dev.off()
    
    grid.arrange(gb,ncol=2)
    plot(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    
    
  }
  
  if(print_to_pdf==FALSE && desicion_to_see_plot ==TRUE){
    
    color_cluster_min_size <- as.data.frame(clp$membership)
    cluster_min_sized_names$V1 <- as.numeric(cluster_min_sized_names$V1)
    color_cluster_min_size$try <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        tmp <- df_molten[grep(paste0("^" , x[1] , "$") , df_molten$Cluster) , ]
        tmp[1 , 6]
      }else{
        "waste"
      }
    })
    cluster_min_sized_names$V1 <- as.character(cluster_min_sized_names$V1)
    color_cluster_min_size$`clp$membership` <- as.character(color_cluster_min_size$`clp$membership`)
    color_cluster_min_size$vertex_size <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        3
      }else{
        1
      }
    })
    hc_add_on <- color_cluster_min_size
    par(mfrow=c(1,2))
    frame()
    heatmap<-Heatmap(dat_new_dendo, 
                     cluster_rows = T,
                     cluster_columns = T,
                     color = new_col)
    gb = grid.grabExpr(draw(heatmap))#color auch hier ersetzt      pdf(name_for_netowrk_pdf , width = 8 , height = 8)
    

    grid.arrange(gb,ncol=2)
    
    Sys.sleep(5)
    V(g)$size<-color_cluster_min_size$vertex_size
    color_cluster_min_size$vertex_color <- color_cluster_min_size$try
    color_cluster_min_size$vertex_color <- apply(color_cluster_min_size,1, function(x){
      
      if(x[4]=="waste"){"white"}else{x[4]}
      
    })
    #color_cluster_min_size[color_cluster_min_size$try%in%c("waste"),]<-c("grey")
    plot(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    
  }
  
  if(print_to_pdf==TRUE && desicion_to_see_plot ==FALSE){
    
    color_cluster_min_size <- as.data.frame(clp$membership)
    cluster_min_sized_names$V1 <- as.numeric(cluster_min_sized_names$V1)
    color_cluster_min_size$try <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        tmp <- df_molten[grep(paste0("^" , x[1] , "$") , df_molten$Cluster) , ]
        tmp[1 , 6]
      }else{
        "waste"
      }
    })
    cluster_min_sized_names$V1 <- as.character(cluster_min_sized_names$V1)
    color_cluster_min_size$`clp$membership` <- as.character(color_cluster_min_size$`clp$membership`)
    color_cluster_min_size$vertex_size <- apply(color_cluster_min_size , 1 , function(x) {
      if(x[1] %in% cluster_min_sized_names$V1 & !(x[1]==c("0"))) {
        3
      }else{
        1
      }
    })
    hc_add_on <- color_cluster_min_size
    par(mfrow=c(1,2))
    frame()
    heatmap<-Heatmap(dat_new_dendo, 
                     cluster_rows = T,
                     cluster_columns = T,
                     color = new_col)
    gb = grid.grabExpr(draw(heatmap))#color auch hier ersetzt      pdf(name_for_netowrk_pdf , width = 8 , height = 8)
    
    pdf(paste0(name_for_pdf,".pdf") , width = 8 , height = 8)
    
    grid.arrange(gb,ncol=2)
    
    Sys.sleep(5)
    V(g)$size<-color_cluster_min_size$vertex_size
    color_cluster_min_size$vertex_color <- color_cluster_min_size$try
    color_cluster_min_size$vertex_color <- apply(color_cluster_min_size,1, function(x){
      
      if(x[4]=="waste"){"white"}else{x[4]}
      
    })
    #color_cluster_min_size[color_cluster_min_size$try%in%c("waste"),]<-c("grey")
    plot(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    
    #ppp<-ggnet2(g , vertex.color = color_cluster_min_size$vertex_color , layout = layout_for_network , vertex.frame.color = "lightgrey" , vertex.label = NA) #, vertex.size = as.numeric(color_cluster_min_size$vertex_size))
    dev.off()
    
  }
  
  color_cluster_min_size$Gene<-clp$names
  RT<-list()
  RT[["color_cluster_min_size"]]<-color_cluster_min_size
  RT[["summary"]]<-summary
  RT[["heatmap_df"]]<-dat_new_dendo
  return(RT)
}


cluster_to_excel<-function(data=cluster_information,sample_cluster = T){
  
  newpath <- paste0(cutoff_wd,"/heatmaps")
  dir.create(file.path(newpath))
  setwd(file.path(newpath))
  
  #triv_gene <- getBM(attributes = c("entrezgene","wikigene_description"), mart = mouse)
  Excel <- xlsx::createWorkbook(type="xlsx")
  for(cluster in unique(data$try)){
    dev.off()
    print(paste0("Start: ",cluster))
    sheet <- xlsx::createSheet(Excel, sheetName = cluster)
    genes<-data[data$try==cluster,"Gene"]
    print(paste0("N genes: ",length(genes)))
  
    GFC_cluster_genes<-GFC_all_genes[genes,]
    testi<-mygene::queryMany(genes,scopes = "symbol",species=organism,returnall=TRUE,return.as = "DataFrame")
    duplicates<-gsub("X","",colnames(testi$duplicates))
    missing<-testi$missing
    testi<-as.data.frame(testi$response)
    testi<-testi[!(grepl("ENS",testi$X_id)),]
    colnames(testi)[1]<-"Gene"
    testi$X_id<-NULL
    testi$X_score<-NULL
    output_table<-merge(testi,GFC_cluster_genes,by="Gene")
    expression_data_filter<-original_data[genes,]
    expression_data <- data.frame(t(expression_data_filter),check.names = FALSE)
    expression_data$Mean_group<-info_Dataset[,"merged"]
    expression_data <- expression_data[ , c(ncol(expression_data) , 1:(ncol(expression_data)-1))]
    mean_expression_per_condition <- setNames(data.frame(t(expression_data[ , -1])) , expression_data[,1])
    mean_expression_per_condition <- t(apply(mean_expression_per_condition , 1 , function(x) tapply(x , colnames(mean_expression_per_condition) , mean)))
    mean_expression_per_condition<-as.data.frame(mean_expression_per_condition)
    colnames(mean_expression_per_condition)<-paste0(colnames(mean_expression_per_condition),"-mean")
    mean_expression_per_condition$Gene<-rownames(mean_expression_per_condition) 
    expression_data_filter$Gene<-rownames(expression_data_filter)
    tmp<-merge(mean_expression_per_condition,expression_data_filter,by="Gene")
    output_table<-merge(output_table,tmp,by="Gene")
    print(paste0("excel name: ",cluster))
    
    xlsx.addTable(Excel, sheet, output_table)
    expression_data_filter$Gene<-NULL
    pdf(paste0(cluster,"_Heatmap.pdf"))
    print(paste0("heatmap name: ",cluster))
    
    pheatmap::pheatmap(expression_data_filter, 
                       cluster_row = T,
                       cluster_cols = sample_cluster,
                       color = col.pal,
                       scale = c("row"),
                       annotation_col = info_Dataset,
                       main = c("Gene Expression Data - whole Dataset "))
    #dev.off() 
  }
  xlsx::saveWorkbook(Excel, "Cluster_information.xlsx")
  dev.off()
  setwd(file.path(cutoff_wd))
}


##GFC calculation
GFC_calculation<-function(normdata=original_data,
                          group=c("merged"),
                          data_in_log=FALSE,
                          range_GFC=2.5){
  # summary[[c("GFC_calculation")]]<-data.frame(group,data_in_log,range_GFC)
  summary<-allargs()
  norm_data_anno <- data.frame(t(normdata),check.names = FALSE)
  
  if(!(group %in% colnames(info_Dataset))){
    group<-colnames(info_Dataset)[1]
    print(paste0("'",colnames(info_Dataset)[1],"'"," was chosen as group variable for GFC - you can specifiy this in the function"))
  }
  
  if(length(group)>1){
    tmp_merged<-data.frame(1:nrow(info_Dataset))
    for(column in group){
      tmp_merged[,column]<-info_Dataset[,column]
    }
    rownames(tmp)<-rownames(info_Dataset)
    tmp_merged$merged<-apply(tmp_merged,1,function(x){
      paste0(x,collapse = "-")
    })
    norm_data_anno$GFC_group<- tmp_merged$merged ##eventuell mergen um fehler vorzubeugen ??
  }else{
    norm_data_anno$GFC_group<- info_Dataset[,group]
  }
  norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
  trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])
  
  if(data_in_log==TRUE){
    antilog <- function(lx , base) {
      
      lbx <- lx/log(exp(1) , base = base)
      result <- exp(lbx)
      result
    }
    trans_norm <- antilog(trans_norm , 2)
  }
  trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , mean)))
  trans_norm <- cbind(trans_norm , rowMeans(trans_norm))
  colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
  trans_norm_temp <- trans_norm
  
  for(col_elements in colnames(trans_norm_temp[ , 1:ncol(trans_norm_temp)-1])) {
    GFC <- foldchange(trans_norm_temp[ , col_elements] , trans_norm_temp[ , ncol(trans_norm_temp)])
    trans_norm <- cbind(trans_norm , GFC)
    colnames(trans_norm)[ncol(trans_norm)] <- paste("GFC_" , col_elements , sep = "") 
  }
  
  trans_norm <- round(trans_norm , digits = 3)
  trans_norm <- as.data.frame(trans_norm)
  GFC_all_genes <- trans_norm[ , grepl("GFC" , colnames(trans_norm))] ##reutnr!
  
  limitsup <- range_GFC
  limitsdown <- (-range_GFC)
  for(iii in 1:(ncol(GFC_all_genes))) {
    GFC_all_genes[ , iii] <- apply(GFC_all_genes[iii] , 1 , function(x) ifelse(x > (limitsup) , limitsup , x)) 
    GFC_all_genes[ , iii] <- apply(GFC_all_genes[iii] , 1 , function(x) ifelse(x < (limitsdown) , (limitsdown) , x)) 
  }
  GFC_all_genes$Gene <- rownames(GFC_all_genes)
  RT<-list()
  RT[["GFC_all_genes"]]<-GFC_all_genes
  RT[["summary"]]<-summary
  
 
  return(RT)
}

##toCytoscape

toCytoscape_all<-function(cluster_information=cluster_information){
  cluster_information_df<-cluster_information
  
  if(organism=="human"){
    orgi = "Human"
  }else{
    orgi="Mouse"
  }
  
  if(checkCytoscapeVersion()[2]==c("3.5.1")){
    port.number = 1234
    base.url = paste('http://localhost:' , port.number , '/v1' , sep = "")
  }else{
    print("Skript is tested for Version 3.5.1 , may work with higher versions")
    port.number = 1234
    base.url = paste('http://localhost:' , port.number , '/v1' , sep = "")
  }
  cluster_information_df$vertex_size<-cluster_information_df$vertex_size*8
  cluster_information_df$TF_list <- apply(cluster_information_df , 1 , function(x) {
    if(x["Gene"] %in% TF_list[,orgi]) {                                       
      tmp <- TF_list[grepl(paste0("^",x["Gene"],"$") , TF_list[,orgi]) , ]
      as.character(tmp[1 , 3])
    }else{
      c(" ")
    }})
  cluster_information_df$TF_Info_border <- apply(cluster_information_df , 1 , function(x) {
    if(x["Gene"] %in% TF_list[,orgi]) {                                      
      c("red")
    }else{
      c("grey")
    }})
  cluster_information_df$TF_label<- apply(cluster_information_df,1,function(x){
    if(x["Gene"]%in% TF_list[,orgi]){
      x["Gene"]
    }else{
      c("")
    }
  })
  
  
  cluster_information_df<-merge(cluster_information_df,GFC_all_genes,by="Gene")
  colnames_with_GFC<-colnames(cluster_information_df)[grepl("GFC",colnames(cluster_information_df))]
  for(col in colnames_with_GFC){
    plot <- ggplot(cluster_information_df , aes(x = rownames(cluster_information_df) , y = cluster_information_df[,col]))+
      geom_point(aes(color = as.numeric(cluster_information_df[,col])))+
      scale_colour_gradientn(colours = c("red","white", "blue"),
                             limits = c(max(cluster_information_df[,col]) , min(cluster_information_df[,col])))
    gg <- ggplot_build(plot)
    colname<- paste0(col,"_Color")
    cluster_information_df[ , colname] <- gg$data[[1]]["colour"] 
  }
  rownames(cluster_information_df)<-cluster_information_df$Gene
  cluster_information_df<-cluster_information_df[V(igraph_object)$name,]
  cluster_information_df$Gene<-NULL
  layout<-as.data.frame(layout)
  colnames(layout)<-c("X","Y")
  cluster_information_df<-cbind(cluster_information_df,layout)
  colnames(cluster_information_df)<-gsub("-","_",colnames(cluster_information_df))
  mynodes <- data.frame(id = as.character(V(igraph_object)$name) , cluster_information_df , stringsAsFactors = FALSE)
  myedges <- data.frame(source = as.character(get.edgelist(igraph_object)[ , 1]) , target = as.character(get.edgelist(igraph_object)[ , 2]),color=c("lightgrey") , stringsAsFactors = FALSE)
  network.name = c("myNetwork")
  collection.name = "myCollection"
  network.suid <- createNetwork(mynodes , myedges , network.name , collection.name)
  
  #ig <- graph_from_data_frame(myedges, directed=F, vertices=mynodes)
  #network.suid <- createNetworkFromIgraph(ig,"myIgraph")
  
  ###delete all styles
  #listStyles(base.url)
  style_names<-listStyles(base.url)
  for(style in style_names){
    style.url = paste(base.url, "styles", sep="/")
    style.delete.url = paste(style.url,style, sep="/")
    DELETE(url=style.delete.url)
  }
  
  style.name = "default_Network"
  defaults <- list(NODE_SHAPE = "Ellipse" , NODE_SIZE = 20 , NODE_FILL_COLOR = "lightgreen")
  nodePositionX <- mapVisualProperty("Node X Location" , "X" , "p")
  nodePositionY <- mapVisualProperty("Node Y Location" , "Y" , "p")
  edge_stroke_color <- mapVisualProperty("Edge Color","color","p")
  createStyle(style.name , defaults , list(nodePositionX , nodePositionY,edge_stroke_color))
  
  colnaming<-colnames(cluster_information_df)[grepl("Color",colnames(cluster_information_df))]
  colnaming<-gsub(" ",".",colnaming)
  
  
  for(styles in colnaming) {
    style.name = styles
    defaults <- list(NODE_SHAPE = "Ellipse" , NODE_SIZE = 20 , EDGE_TRANSPARENCY = 80)
    nodeColor <- mapVisualProperty("Node Fill Color" , styles , "p")
    createStyle(style.name , defaults , list( nodeColor,edge_stroke_color))
  }
  
  style.name = "TF_labelled_names"
  defaults <- list(NODE_SHAPE = "Ellipse" , NODE_SIZE = 5 , EDGE_TRANSPARENCY = 20 , NODE_FILL_COLOR = "white" , NODE_LABEL_FONT_SIZE = 10 , NODE_TRANSPARENCY = 255)
  TF_labelled <- mapVisualProperty("Node Label" , "TF_label" , "p")
  TF_fill_color<- mapVisualProperty("Node Fill Color" , "TF_list" , "d",c("TF","Co_factor","Chromatin_remodeller"),c("#99CCFF","#ed9ce7","#FF7777"))
  TF_size<- mapVisualProperty("Node Size" , "vertex_size" , "p")
  TF_Info_border <- mapVisualProperty("Node Border Paint" , "TF_Info_border" , "p")
  
  createStyle(style.name , defaults , list(TF_labelled , TF_fill_color,
                                           TF_size,TF_Info_border,edge_stroke_color))
  
  style.name="cluster_color"
  defaults<-list(NODE_SHAPE="Ellipse",
                 NODE_SIZE=20,
                 EDGE_TRANSPARENCY=20,
                 NODE_FILL_COLOR="white",
                 NODE_TRBEL_FONT_SIZE=10,ANSPARENCY=255)
  cluster_color<-mapVisualProperty("Node Fill Color","try","p")
  createStyle(style.name, defaults, list(cluster_color,edge_stroke_color))
}

LayoutfromCytoscape<-function(){
  port.number = 1234
  base.url = paste('http://localhost:' , port.number , '/v1' , sep = "")
  #setwd(cutoff_wd)
  out_put<-exportNetwork("cytoscape_object", 'CYJS', base.url = base.url)
  out_put_split<-strsplit(out_put,":")[[1]]
  out_put_sd_folder<-gsub("[\\]","/",out_put_split[3])
  out_put_url<-paste0(strsplit(out_put_split[2]," ")[[1]][2],":",out_put_sd_folder)
  my.JSON <- rjson::fromJSON(file=out_put_url)
  file.remove(out_put_url)
  position_list<-list()
  for(list_platz in 1:length(my.JSON$elements$nodes)){
    my.JSON$elements$nodes[[list_platz]]$position
    position_list[[list_platz]]<-data.frame(Gene=my.JSON$elements$nodes[[list_platz]]$data$name,x=as.numeric(my.JSON$elements$nodes[[list_platz]]$position$'x'),y=as.numeric(my.JSON$elements$nodes[[list_platz]]$position$'y'))
  }
  coordinates_cy<-list.rbind(position_list)
  rownames(coordinates_cy)<-coordinates_cy$Gene
  coordinates_cy<-coordinates_cy[cluster_information$Gene,2:3]
  rownames(coordinates_cy)<-1:nrow(coordinates_cy)
  is.matrix(coordinates_cy)
  rownames(coordinates_cy)<-NULL
  colnames(coordinates_cy)<-NULL
  return(coordinates_cy)
}



##
test_layout<-function(layouts_to_test=c("layout_with_fr","layout_with_kk","layout_with_lgl"),min_nodes_number_for_network=10, pdf = T){
  # summary[[c("test_layout")]]<-data.frame(layouts_to_test,min_nodes_number_for_network)
  summary<-allargs()
  data<-correlation_df
  cutted_data<-data[data$correlation>chosen_cutoff,]
  igraph_object<-igraph::graph_from_data_frame(cutted_data,directed=FALSE)
  layouts_on_list_coord<-list()
  igraph_list<-list()
  for(layouts_test in layouts_to_test ){
    if(igraph::count_components(igraph_object)==1){
      igraph_list[[1]]<-igraph_object
      layout_df<-BiocGenerics::do.call(layouts_test, igraph_list)
      layouts_on_list_coord[[layouts_test]]<-layout_df
      g<-igraph_object
    }else{
      components<-1:igraph::components(igraph_object)[[3]]
      components_numbers<-components[igraph::components(igraph_object)[[2]]>min_nodes_number_for_network]
      V(igraph_object)$comp<-igraph::components(igraph_object)[[1]]
      list_of_components_graphs<-list()
      
      for(component in components_numbers){
        list_of_components_graphs[[as.character(component)]] <- igraph::induced_subgraph(igraph_object , V(igraph_object)$comp == component)
      }
      layouts_on_list <- lapply(list_of_components_graphs , layouts_test)
      layout_df <- merge_coords(list_of_components_graphs , layouts_on_list)
      g <- disjoint_union(list_of_components_graphs)
      layouts_on_list_coord[[layouts_test]]<-layout_df
    }
  }
  
  if(pdf){
  pdf("layouts.pdf")
  par(mfrow=c(2,2))
  for(listenplatz in names(layouts_on_list_coord)){
    plot(g,layout=layouts_on_list_coord[[listenplatz]],
         vertex.size = 3 ,
         vertex.frame.width = 3/100 ,
         vertex.label = NA , 
         vertex.frame.color = "lightgrey",
         main=listenplatz)
  }
  dev.off()
  }
  else{
    par(mfrow=c(2,2))
    for(listenplatz in names(layouts_on_list_coord)){
      plot(g,layout=layouts_on_list_coord[[listenplatz]],
           vertex.size = 3 ,
           vertex.frame.width = 3/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey",
           main=listenplatz)
    }
  }
  RT<-list()
  RT[["layouts_on_list_coord"]]<-layouts_on_list_coord
  RT[["summary"]]<-summary
  return(RT)
}
###
search_for_good_cutoff<-function(data=correlation_df,min_nodes_number_for_network=10,
                                 show_network=TRUE,layout=layout_with_kk,cluster_algo= cluster_louvain,
                                 max_cluster_count_per_gene=8,
                                 min_cluster_size=10,
                                 abberation=0.1,no_cutoffs_tested=4 ,cutoffs_to_test=c("0")){
  #layout_probem /cluster problem
  
  # summary[[c("search_for_good_cutoff")]]<-data.frame(min_nodes_number_for_network,show_network,max_cluster_count_per_gene,
  #                                           min_cluster_size,abberation,no_cutoffs_tested,cutoffs_to_test)
  
  if(abberation==0){
    range<-as.numeric(cutoffs_to_test)
  }else{
    range<-seq(from = chosen_cutoff*(1-0.1) , to = chosen_cutoff*(1+0.1) , length.out = no_cutoffs_tested)
  }
  i<-0
  heatmap_list<-list()
  for (cutoff in range){
    i<-i+1
    cutted_data<-data[data$correlation>cutoff,]
    igraph_object<-graph_from_data_frame(cutted_data,directed=FALSE)
    
    if(count_components(igraph_object)==1){
      layout_df<-layout(igraph_object)
      g<-igraph_object
      
      
    }else{
      components<-1:igraph::components(igraph_object)[[3]]
      components_numbers<-components[igraph::components(igraph_object)[[2]]>min_nodes_number_for_network]
      V(igraph_object)$comp<-igraph::components(igraph_object)[[1]]
      list_of_components_graphs<-list()
      for(component in components_numbers){
        list_of_components_graphs[[as.character(component)]] <- igraph::induced_subgraph(igraph_object , V(igraph_object)$comp == component)
      }
      layouts_on_list <- lapply(list_of_components_graphs , layout)
      layout_df <- merge_coords(list_of_components_graphs , layouts_on_list)
      g <- igraph::disjoint_union(list_of_components_graphs)
      
    }
    clp<-cluster_algo(g)
    gene_which_cluster<-list()
    gene_which_cluster[[1]]<-as.data.frame(clp$membership)
    gene_which_cluster<-as.data.frame(gene_which_cluster[[1]])
    colnames(gene_which_cluster)<-c("cluster_max")
    gene_which_cluster$cluster_max <- as.numeric(gene_which_cluster$cluster_max)
    for_printing<-as.data.frame(gene_which_cluster$cluster_max,check.names=FALSE)
    p<-as.data.frame(table(for_printing[,1] == 0))
    print(paste0(p[grepl("TRUE",p$Var1),2]," Genes are advised to more than ",max_cluster_count_per_gene," Clusters, and will be painted white and left out in upcoming plots and analysis"))
    clp$membership <- gene_which_cluster$cluster_max  
    gene_which_cluster$cluster_max <- as.numeric(gene_which_cluster$cluster_max)
    
    clp$membership <- gene_which_cluster$cluster_max 
    
    big_cluster <- list()
    cluster_min_sized_names <- list()
    big_cluster_names <- paste0("Cluster" , unique(clp$membership))
    
    for(big_clust in 1:length(unique(clp$membership))){
      if (length(clp[[big_clust]]) >= min_cluster_size) {                                            # here adjust the minimal size of a cluster that you still would like to include
        big_cluster[[big_cluster_names[big_clust]]] <- as.data.frame(clp[[big_clust]])
        cluster_min_sized_names[[big_clust]] <- big_clust
      }
    }
    big_cluster_all <- as.data.frame(list.rbind(big_cluster))
    cluster_min_sized_names <- as.data.frame(list.rbind(cluster_min_sized_names))
    
    list_of_means <- list()
    list_of_GFC_per_cluster<-list()
    names_of_GC_per_cluster <- paste0("Cluster" , 1:length(unique(clp$membership)))
    for(cluster in 1:length(unique(clp$membership))) {
      list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] <- GFC_all_genes[grepl(paste0(paste0("^" , clp[[cluster]] , "$") , collapse = "|") , GFC_all_genes$Gene) , ]
      list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]][ , ncol(GFC_all_genes)] <- NULL
      list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] <- rbind(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] , apply(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] , 2 , mean))
      list_of_means[[cluster]] <- apply(list_of_GFC_per_cluster[[names_of_GC_per_cluster[cluster]]] , 2 , mean)
    }
    list_of_means <- as.data.frame(list.rbind(list_of_means))
    
    list_of_means_new <- list_of_means[cluster_min_sized_names$V1 , ]
    row.order <- hclust(dist(list_of_means_new))$order
    col.order <- hclust(dist(t(list_of_means_new)))$order
    dat_new <- list_of_means_new[row.order , col.order]
    
    df_molten <- reshape2::melt(as.matrix(dat_new))
    names(df_molten) <- c("Cluster" , "Tissue" , "GFC")
    df_molten$Cluster <- as.numeric(df_molten$Cluster)
    df_molten$Gen_No_cluster <- apply(df_molten , 1 , function(x){
      x <- as.numeric(x[1])
      length(clp[[x]])
    })
    df_molten$Gen_No_cluster <- as.numeric(df_molten$Gen_No_cluster)
    df_molten$NEW_1 <- paste0(df_molten$Cluster , "[" , df_molten$Gen_No_cluster , "]")
    df_molten$Colors <- labels2colors(df_molten$Cluster , zeroIsGrey = TRUE)
    df_molten$NEW_2 <- paste0(df_molten$Colors , "[" , df_molten$Gen_No_cluster , "]")
    df_molten$Vertex_size <- as.numeric(c("3"))
    
    dat_new_dendo<-dat_new
    rownames(dat_new_dendo)<-unique(df_molten$NEW_2)
    col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
    
    if(show_network==TRUE){
      color_cluster_min_size <- as.data.frame(clp$membership)
      cluster_min_sized_names$V1 <- as.numeric(cluster_min_sized_names$V1)
      color_cluster_min_size$try <- apply(color_cluster_min_size , 1 , function(x) {
        if(x[1] %in% cluster_min_sized_names$V1) {
          tmp <- df_molten[grep(paste0("^" , x[1] , "$") , df_molten$Cluster) , ]
          tmp[1 , 6]
        }else{
          "white"
        }
      })
      cluster_min_sized_names$V1 <- as.character(cluster_min_sized_names$V1)
      color_cluster_min_size$`clp$membership` <- as.character(color_cluster_min_size$`clp$membership`)
      color_cluster_min_size$vertex_size <- apply(color_cluster_min_size , 1 , function(x) {
        if(x[1] %in% cluster_min_sized_names$V1) {
          3
        }else{
          1.5
        }
      })
      hc_add_on <- color_cluster_min_size
      
      heatmap<-pheatmap::pheatmap(dat_new_dendo, 
                                  cluster_row = T,
                                  cluster_cols = T,
                                  color = col.pal,
                                  title=paste0("Cutoff: ",cutoff),
                                  silent = TRUE)
      
      par(mfrow=c(1,2))
      frame()
      grid.arrange(heatmap[[4]],ncol=2,nrow=1)
      plot(g , vertex.color = color_cluster_min_size$try , layout = layout_df , vertex.frame.color = "lightgrey" ,
           vertex.label = NA , vertex.size = color_cluster_min_size$vertex_size,main=paste0("Cutoff: ",cutoff))
      
      
    }else{
      
      
      
      heatmap<-pheatmap::pheatmap(dat_new_dendo, 
                                  cluster_row = T,
                                  cluster_cols = T,
                                  color = col.pal,
                                  main=paste0("Cutoff: ",cutoff),
                                  silent = TRUE)
      heatmap_list[[i]]<-heatmap[[4]]
      
    }
    
  }
  if(show_network==FALSE){
    grid<-do.call(grid.arrange,heatmap_list)
    grid
  }
}

##plot GFC
##plot GFC Networks
plot_GFC_networks<-function(igraph_object=igraph_object,print_to_pdf=FALSE,print_edges_png_nodes_pdf=c("each","one","none"),GFC_all_genes=GFC_all_genes){
  # summary[["plot_GFC_networks"]]<-data.frame(print_to_pdf)
  
  GFC_all_genes_color<-GFC_all_genes
  
  layout_df<-cbind(layout,Gene=V(igraph_object)$name)
  GFC_all_genes_color<-merge(GFC_all_genes_color,layout_df,by="Gene")
  
  if(is.null(layout)){
    print("Make sure your layout varaible has been definied!")
  }
  
  GFC_all_genes_color$V1<-as.numeric(as.character(GFC_all_genes_color$V1))
  GFC_all_genes_color$V2<-as.numeric(as.character(GFC_all_genes_color$V2))
  
  for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
    library(ggplot2)
    plot <- ggplot2::ggplot(GFC_all_genes_color,aes(x = V1 , y = V2))+
      geom_point(aes(colour = as.numeric(GFC_all_genes_color[,col])))+
      scale_colour_gradientn(colours = c("blue","white", "red"),
                             limits = c(min(GFC_all_genes_color[,col]) , max(GFC_all_genes_color[,col])))
    gg <- ggplot2::ggplot_build(plot)
    #colname<- paste0(col,"_Color")
    GFC_all_genes_color[ , col] <- gg$data[[1]]["colour"] 
  }
  rownames(GFC_all_genes_color)<-GFC_all_genes_color$Gene
  GFC_all_genes_color<-GFC_all_genes_color[V(igraph_object)$name,]
  if(print_to_pdf==FALSE){
    par(mfrow = c(2 , ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2))))
    for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
      igraph::plot.igraph(igraph_object,layout=layout,
           vertex.size = 3 ,
           vertex.frame.width = 0.01/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey",
           vertex.color=GFC_all_genes_color[,col],
           main=col)
    }
  }else{
    mainDir <- originalwd
    subDir <- paste0(chosen_cutoff,"/GFC")
    dir.create(file.path(mainDir , subDir))
    setwd(file.path(mainDir , subDir))
    cutoff_wd<-getwd()
    pdf("GFC_colord_network.pdf" , width = 8 , height = 8)
    par(mfrow = c(ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2)),2))
    for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
      igraph::plot.igraph(igraph_object,layout=layout,
           vertex.size = 3 ,
           vertex.frame.width = 0.01/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey",
           vertex.color=GFC_all_genes_color[,col],
           main=col)
    }
    dev.off()
    
  }
  if(print_edges_png_nodes_pdf==c("one")){
    mainDir <- originalwd
    subDir <- paste0(chosen_cutoff,"/GFC")
    dir.create(file.path(mainDir , subDir))
    setwd(file.path(mainDir , subDir))
    cutoff_wd<-getwd()
    pdf("GFC_colord_network_nodes.pdf" , width = 8 , height = 8)
    par(mfrow = c(ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2)),2))
    for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
      igraph::plot.igraph(igraph_object,layout=layout,
           vertex.size = 3 ,
           vertex.frame.width = 0.01/100 ,
           vertex.label = NA , 
           vertex.frame.color = "lightgrey",
           vertex.color=GFC_all_genes_color[,col],
           main=col,
           edge.width=NA)
      
    }
    dev.off()
    
    png("GFC_colord_network_edges.png" , width = 8 , height = 8,res=680,units = "in")
    par(mfrow = c(ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2)),2))
    for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
      igraph::plot.igraph(igraph_object,layout=layout,
           vertex.size = 3 ,
           vertex.frame.width = NA ,
           vertex.label = NA , 
           vertex.frame.color = NA,
           vertex.color=NA,
           main=NA,
           edge.width=0.1)
      
    }
    dev.off()
    
  }else{
    if(print_edges_png_nodes_pdf==c("each")){
      mainDir <- originalwd
      subDir <- paste0(chosen_cutoff,"/GFC")
      dir.create(file.path(mainDir , subDir))
      setwd(file.path(mainDir , subDir))
      cutoff_wd<-getwd()
      
      for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
        pdf(paste0(col,"_nodes.pdf") , width = 8 , height = 8)
        #par(mfrow = c(2 , ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2))))
        igraph::plot.igraph(igraph_object,layout=layout,
             vertex.size = 3 ,
             vertex.frame.width = 0.01/100 ,
             vertex.label = NA , 
             vertex.frame.color = "lightgrey",
             vertex.color=GFC_all_genes_color[,col],
             main=col,
             edge.width=NA)
        dev.off()
      }
      
      for(col in colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))]){
        png(paste0(col,"_edges.png") , units ="in", width = 8 , height = 8,res=680)
        #par(mfrow = c(2 , ceiling((length(colnames(GFC_all_genes_color)[grepl("GFC", colnames(GFC_all_genes_color))])/2))))
        igraph::plot.igraph(igraph_object,layout=layout,
             vertex.size = 3 ,
             vertex.frame.width = NA ,
             vertex.label = NA , 
             vertex.frame.color = NA,
             vertex.color=NA,
             main=NA,
             edge.width=0.1)
        dev.off()
      }
      
      
    }
  
  }
}



#####

plotFunction<-function(){
  p<- ggplot(TFs_matrix_all_plot.melt.mean,aes(x=variable,y=mean,fill=merged)) +
    geom_bar(stat="identity",position = "dodge") +
    facet_wrap(~variable,scales = "free")
  print(p)
}

#clusterprofiler
clusterprofiler_autoCena<-function(cluster_to_check=c("all",cluster_color),group=c("merged")){
  summary<-allargs()
  cluster_available<-unique(cluster_information$try)
  df<-data.frame()
  no<-0
  biggest_cluster <- cluster_information[grepl(paste0("^" , names(which.max(table(cluster_information$try))) , "$") , cluster_information$try) , ]
  df <- data.frame(longest_cluster = as.character(cluster_information[cluster_information$try == unique(biggest_cluster$try) , 2]))
  
  if(cluster_to_check==c("all")){
    for(Thomas in cluster_available) {
      no <- 1 + no
      x <- as.character(Thomas)
      new.col <- as.character(cluster_information[cluster_information$try == x , 5])
      df[ , no] <- c(new.col , rep(NA , nrow(df)-length(new.col)))
      colnames(df)[no] <- Thomas
    }
  }else{
    for(Thomas in cluster_to_check){
      no <- 1 + no
      x <- as.character(Thomas)
      new.col <- as.character(cluster_information[cluster_information$try == x , 5])
      df[ , no] <- c(new.col , rep(NA , nrow(df)-length(new.col)))
      colnames(df)[no] <- Thomas
    }
  }
  
  ###cluster Profiler
  normdata<- as.data.frame(t(Dataset_1))
  normdata$merged<-info_Dataset[,group]
  universe <- rownames(Dataset_1)
  
  cluster_genes_original<-data.frame(Gene=df)
  returnerlist<-list()
  
  if(!(base::exists(c("cutoff_wd")))){
    cutoff_wd<-file.path(originalwd,chosen_cutoff)
  }
  
  x <- cutoff_wd
  dir.create(x)
  setwd(x)
  plotPath = file.path(x, "clusterProfiler");
  dir.create(file.path(x, "clusterProfiler"), showWarnings = FALSE)
  
  if(organism==c("Mouse")){
    orga_type= c("mouse")
  }else{
    orga_type=c("human")
  }
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("biomaRt")
  #library(biomaRt)
  
  #human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  if(orga_type == "mouse"){
    universe_mouse_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = universe, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    universe_mouse_human <- universe_mouse_human[,2]
    universe_Entrez_mouse = clusterProfiler:: bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    universe_Entrez_mouse = unlist(universe_Entrez_mouse[2],use.names = FALSE)
    universe_Entrez_mouse_human = clusterProfiler::bitr(universe_mouse_human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez_mouse_human = unlist(universe_Entrez_mouse_human[2],use.names = FALSE)
  }else{
    universe_Entrez = clusterProfiler:: bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez = unlist(universe_Entrez[2],use.names = FALSE)
  }
  c1_hallmark_genes <- gmtfile
  list_of_entrez <- list()
  #######
  for(id_it in 1:ncol(cluster_genes_original)){
    

    list_of_genes <- list(cluster_genes_original[!is.na(cluster_genes_original[,id_it]),id_it])

    print(paste0("Prediction for cluster: ",colnames(cluster_genes_original)[id_it]," [",length(!is.na(unlist(list_of_genes))),"]"))
    
    
    #if mouse, mouse symbols are translated to human, because some gene sets works only with human
    if(orga_type == "mouse"){
      universe_orignal <- universe
      genes_cluster = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unlist(list_of_genes), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      cluster_genes <- genes_cluster[,2]
      list_of_genes_mouse_human <- list(cluster_genes)
      entrez_de = clusterProfiler :: bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
      module_entrez_mouse <- unlist(entrez_de[2],use.names = FALSE)
      entrez_de = clusterProfiler::bitr(unlist(list_of_genes_mouse_human), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez_mouse_human <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez_mouse
      
    }else{
      entrez_de = clusterProfiler::bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez
    }
    
    wb <- xlsx::createWorkbook(type="xlsx")

    #wb <- createWorkbook()
    
    
    pdf(paste("clusterProfiler/ClusterProfiler_",colnames(cluster_genes_original)[id_it],".pdf",sep=""), onefile=FALSE,  width = 15, height = 15)
    #pdf(paste("ClusterProfiler_",colnames(cluster_genes_original)[id_it],".pdf",sep=""), onefile=FALSE,  width = 15, height = 15)
    font_size <- 8
    
    
    #Create figure window and layout
    plot.new()
    grid::grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 2)))
    
    #Hallmark gene sets
    #start
    #egmt <- enricher(module_entrez, TERM2GENE=c1_hallmark, universe = universe_genes, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    
    if(orga_type == "mouse"){
      egmt <- clusterProfiler::enricher(unlist(list_of_genes_mouse_human), TERM2GENE=c1_hallmark_genes, universe = universe_mouse_human, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    }else{
      egmt <- clusterProfiler::enricher(unlist(list_of_genes), TERM2GENE=c1_hallmark_genes, universe = universe, pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1.0)
    }
    
    if(!is.null(egmt)){
      hallmark_plot <- dotplot(egmt, font.size = font_size, title = "Hallmark enrichment")
      
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
      print(hallmark_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Hallmark")
      #writeData(wb, 1, egmt@result)
      
      sheet <- xlsx::createSheet(wb, sheetName = "Hallmark")
      xlsx.addTable(wb, sheet, egmt@result)
    }
    #ende
    
    #KEGG
    #start
    #library(clusterProfiler)
    if(orga_type == "mouse"){
      KEGG <- enrichKEGG(gene = module_entrez_mouse, organism = 'mmu', pvalueCutoff=0.05,universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_Kegg<-NULL
      Enriched_Kegg_obj<-NULL
      reg_Mm=AnnotationDbi::select(org.Mm.eg.db,as.character(unlist(list_of_genes)),"ENTREZID","SYMBOL",multiVals="first")
      
      
      if(!is.null(KEGG))
      {
        if(nrow(summary(KEGG))>0)
        {
          df_kk<-as.data.frame(summary(KEGG))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_Kegg<-df_kk
          Enriched_Kegg_obj<-KEGG
        }
        else
        {
          Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_Kegg_obj<-NULL
        }
        
      }
      
    }else{
      KEGG <- enrichKEGG(gene = module_entrez, organism = 'hsa', pvalueCutoff=0.05,universe = universe_Entrez, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_Kegg<-NULL
      Enriched_Kegg_obj<-NULL
      reg_Hs=AnnotationDbi::select(org.Hs.eg.db,as.character(unlist(list_of_genes)),"ENTREZID","SYMBOL",multiVals="first")
      #str(KEGG)
      
      if(!is.null(KEGG))
      {
        if(nrow(summary(KEGG))>0)
        {
          df_kk<-as.data.frame(summary(KEGG))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_Kegg<-df_kk
          Enriched_Kegg_obj<-KEGG
        }
        else
        {
          Enriched_Kegg<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_Kegg_obj<-NULL
        }
        
      }
      
    }
    
    
    if(!is.null(Enriched_Kegg)){
      
      KEGG_plot <- dotplot(KEGG, font.size = font_size, title = "KEGG enrichment")
      pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
      print(KEGG_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "KEGG")
      #writeData(wb, 2, Enriched_Kegg)
      sheet <- xlsx::createSheet(wb, sheetName = "KEGG")
      xlsx.addTable(wb, sheet, Enriched_Kegg)
    }
    #ende
    
    #Reactome
    #start
    if(orga_type == "mouse"){
      ReacTome <- ReactomePA::enrichPathway(gene=module_entrez_mouse, organism = "mouse", pvalueCutoff=0.05, universe = universe_Entrez_mouse, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      
      Enriched_ReacTome<-NULL
      Enriched_ReacTome_obj<-NULL
      
      
      if(!is.null(ReacTome))
      {
        if(nrow(summary(ReacTome))>0)
        {
          df_kk<-as.data.frame(summary(ReacTome))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_ReacTome<-df_kk
          Enriched_ReacTome_obj<-ReacTome@result
        }
        else
        {
          Enriched_ReacTome<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_ReacTome_obj<-NULL
        }
        
      }
      
    }else{
      ReacTome <- ReactomePA::enrichPathway(gene=module_entrez, organism = "human", pvalueCutoff=0.05, universe = universe_Entrez, pAdjustMethod = "none", qvalueCutoff = 1.0)
      
      Enriched_ReacTome<-NULL
      Enriched_ReacTome_obj<-NULL
      
      
      if(!is.null(ReacTome))
      {
        if(nrow(summary(ReacTome))>0)
        {
          df_kk<-as.data.frame(summary(ReacTome))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_ReacTome<-df_kk
          Enriched_ReacTome_obj<-ReacTome@result
        }
        else
        {
          Enriched_ReacTome<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_ReacTome_obj<-NULL
        }
        
      }
    }
    
    
    if(!is.null(Enriched_ReacTome)){
      
      #ReacTome$Description <- factor(ReacTome$Description, levels = rev(unique(ReacTome$Description)))
      
      ReacTome_plot <- dotplot(ReacTome, font.size = font_size, title = "Reactome enrichment")
      #ReacTome_plot <- dotplot(KEGG, font.size = font_size, title = "KEGG enrichment")
      #ReacTome_plot <- dotplot(KEGG, font.size = font_size, title = "Reactome enrichment")
      pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
      print(ReacTome_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Reactome")
      #writeData(wb, 3, Enriched_ReacTome)
      sheet <- xlsx::createSheet(wb, sheetName = "Reactome")
      xlsx.addTable(wb, sheet, Enriched_ReacTome)
    }
    #ende
    
    #enrichDO
    #start
    if(orga_type == "mouse"){
      enrichDO <- enrichDO(gene= module_entrez_mouse_human,
                           ont           = "DO", 
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "none",
                           universe      = universe_Entrez_mouse_human, 
                           qvalueCutoff  = 1.0,
                           readable      = FALSE)
      
      Enriched_enrichDO<-NULL
      Enriched_enrichDO_obj<-NULL
      
      
      if(!is.null(enrichDO))
      {
        if(nrow(summary(enrichDO))>0)
        {
          df_kk<-as.data.frame(summary(enrichDO))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Mm$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Mm$SYMBOL[id], collapse = '/')
          }
          Enriched_enrichDO<-df_kk
          Enriched_enrichDO_obj<-enrichDO@result
        }
        else
        {
          Enriched_enrichDO<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_enrichDO_obj<-NULL
        }
        
      }
      
      
    }else{
      enrichDO <- enrichDO(gene= module_entrez,
                           ont           = "DO", 
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "none",
                           universe      = universe_Entrez, 
                           qvalueCutoff  = 1.0,
                           readable      = FALSE)
      
      Enriched_enrichDO<-NULL
      Enriched_enrichDO_obj<-NULL
      
      
      if(!is.null(enrichDO))
      {
        if(nrow(summary(enrichDO))>0)
        {
          df_kk<-as.data.frame(summary(enrichDO))[1:8]
          for(x in 1:length(df_kk[,8]))
          {
            temp<-strsplit(df_kk[x,8],"/")
            id<-which(reg_Hs$ENTREZID %in% temp[[1]] )
            df_kk[x,8]<-paste(reg_Hs$SYMBOL[id], collapse = '/')
          }
          Enriched_enrichDO<-df_kk
          Enriched_enrichDO_obj<-enrichDO@result
        }
        else
        {
          Enriched_enrichDO<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_enrichDO_obj<-NULL
        }
        
      }
      
    }
    
    
    
    if(!is.null(Enriched_enrichDO)){
      
      enrichDO_plot <- dotplot(enrichDO, font.size = font_size, title = "Disease enrichment")
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 3))
      print(enrichDO_plot, newpage = FALSE)
      popViewport()
      
      #addWorksheet(wb, "Disease")
      #writeData(wb, 4, enrichDO@result)
      sheet <- xlsx::createSheet(wb, sheetName = "Disease")
      xlsx.addTable(wb, sheet, Enriched_enrichDO)
    }
    #ende
    
    #enrichGO_dotplot
    #start
    if(orga_type == "mouse"){
      enrichGO <- clusterProfiler::enrichGO(gene = module_entrez_mouse,
                           universe = universe_Entrez_mouse,
                           OrgDb = org.Mm.eg.db,
                           #keytype = 'ENTREZID',
                           ont = "BP",
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 1,
                           readable      = T)
      
    }else{
      enrichGO <- clusterProfiler::enrichGO(gene = module_entrez,
                           universe = universe_Entrez,
                           OrgDb = org.Hs.eg.db,
                           #keytype = 'ENTREZID',
                           ont = "BP",
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = T)
      
    }
    
    
    if(!is.null(enrichGO)){
      
      enrichGO_plot <- dotplot(enrichGO, font.size = font_size, title = "GO enrichment",showCategory=20)
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 2))
      print(enrichGO_plot, new = FALSE)
      popViewport()
      
      #addWorksheet(wb, "GO")
      #writeData(wb, 5, enrichGO@result)
      sheet <- xlsx::createSheet(wb, sheetName = "GO")
      xlsx.addTable(wb, sheet, enrichGO@result)
    }
    #ende
    
    #enrichGO_enrichment_map
    #start
    
    # pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 2))
    # par(fig = gridFIG(), new = TRUE)
    # enrichMap(enrichGO, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.8, n=10)
    # popViewport()
    #}
    #ende
    # 
    #TF prediction
    plotFunction<-function(){
      p<- ggplot2::ggplot(TFs_matrix_all_plot.melt.mean,aes(x=variable,y=mean,fill=merged)) +
        geom_bar(stat="identity",position = "dodge") +
        facet_wrap(~variable,scales = "free")
      print(p)
    }
    
    
    #start
    if(orga_type == "mouse"){
      
      r <- as.vector(list_of_genes[[1]])
      
      TFtable <- primo(as.vector(list_of_genes[[1]]), inputType = "geneSymbol", org = "Mm")
      #class(list_of_genes)
      
      
      TFs <- str_split(TFtable$overRepresented[,4], "_")
      
      TFs <- do.call(rbind, TFs)[,1]
      
      TFs <- str_split(TFs, "::")
      
      TFs <- do.call(rbind, TFs)
      
      if(grepl(":", TFs)&&!is.null(TFs)){
        
        TFs <- c(TFs[,1],TFs[,2])
        
      }
      
      TFs <-unique(TFs)
      
      TFs <- str_to_title(TFs, locale = "")
      
      
      
    }else{
      
      #TFtable <- primo(unlist(list_of_genes), inputType = "geneSymbol", org = "Hs")
      TFtable <- primo(as.vector(list_of_genes[[1]]), inputType = "geneSymbol", org = "Hs")
      
      
      TFs <- str_split(TFtable$overRepresented[,4], "_")
      
      TFs <- do.call(rbind, TFs)[,1]
      
      TFs <- str_split(TFs, "::")
      
      TFs <- do.call(rbind, TFs)
      
      if(grepl(":", TFs)&&!is.null(TFs)){
        
        TFs <- c(TFs[,1],TFs[,2])
        
      }
      
      TFs <-unique(TFs)
      
      TFs <- str_to_upper(TFs, locale = "")
      
    }  
    
    
    
    
    
    if(!is.null(TFs)){
      
      #addWorksheet(wb, "TFoverrepresented")
      
      #writeData(wb, 6, TFtable$overRepresented)
    
      if(length(TFs) != 0L){
        
        sheet <- xlsx::createSheet(wb, sheetName = "TFoverrepresented")
        
        xlsx.addTable(wb, sheet, TFtable$overRepresented)
        
        TFS_intersect <- intersect(TFs,colnames(normdata))
       
        if(length(TFS_intersect) != 0L){
          
          TFs_matrix_all_plot <- normdata[,c("merged",TFS_intersect)]
          
          TFs_matrix_all_plot.melt <- melt(TFs_matrix_all_plot, id="merged")
          
          # calculate means
          
          require(dplyr)
          library(dplyr)
          TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt %>% dplyr::group_by(merged,variable) %>% dplyr::summarise(mean=mean(value))
          # # get the mean over merged group (90 to 30 samples)
          # TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt %>% group_by(merged,variable)
          # # commas into dots
          # TFs_matrix_all_plot.melt.mean$value <- as.numeric(sub(",", ".", TFs_matrix_all_plot.melt.mean$value, fixed = TRUE))
          # # reorder the results according to group and not TF (!?)
          # TFs_matrix_all_plot.melt.mean <- TFs_matrix_all_plot.melt.mean %>%  dplyr::summarise(mean=mean(value)) 
          
          
          condis <- unique(normdata$merged)
          
          TFs_matrix_all_plot.melt.mean$merged <- factor(TFs_matrix_all_plot.melt.mean$merged,levels = condis)
          
        
          ####plotFunction####
          plotTheThing<-function(){
            
            p<- ggplot(TFs_matrix_all_plot.melt.mean,aes(x=variable,y=mean,fill=merged)) +
              
              geom_bar(stat="identity",position = "dodge") +
              facet_wrap(~variable,scales = "free")
            
            print(p)
            
          }
          
          
          # here you put the plot into the excel sheet
          xlsx.addPlot(wb, sheet, plotTheThing)
          
        }
        
        
        
        TFS_intersect <- intersect(TFs,colnames(normdata))
        
        
        
        if(length(TFS_intersect) != 0L){


          TFs_matrix <- normdata[,c("merged",head(TFS_intersect,4))]


          TFs_matrix.melt <- melt(TFs_matrix, id="merged")        # calculate means

          #require(dplyr)


          TFs_matrix.melt %>% group_by(merged,variable) -> TFs_matrix.melt.mean

          TFs_matrix.melt.mean$value<-as.numeric(sub(",", ".", TFs_matrix.melt.mean$value, fixed = TRUE))

          #TFs_matrix.melt.mean %>% summarise(mean=mean(value)) -> TFs_matrix.melt.mean

          condis <- unique(normdata$merged)

          TFs_matrix.melt.mean$merged <- factor(TFs_matrix.melt.mean$merged,levels = condis)


          # for pdf
          TF_plot <- ggplot(TFs_matrix.melt.mean,aes(x=variable,y=value,fill=merged)) +

            geom_bar(stat="identity",position = "dodge") +

            facet_wrap(~variable,scales = "free")

          # embed in pdf
          pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))

          print(TF_plot, newpage = FALSE)

          popViewport()

        }
        
        
      }
      
    }
    
    xlsx::saveWorkbook(wb, paste("clusterProfiler/ClusterProfiler_",colnames(cluster_genes_original)[id_it],".xlsx",sep=""))
    dev.off()
    
  }
  
  #dev.off()
  
  returnerlist<-list()
  if(!(is.null(egmt))){
  returnerlist[[c("Hallmark")]]<-egmt@result
  }
  
  returnerlist[[c("KEGG")]]<-Enriched_Kegg
  returnerlist[[c("Enriched_ReacTome")]]<-Enriched_ReacTome
  returnerlist[[c("Enriched_enrichDO")]]<-Enriched_enrichDO
  if(!(is.null(enrichGO))){
  returnerlist[[c("enrichGO")]]<- enrichGO@result
  }
  if(!(is.null(TFtable$overRepresented))){
  returnerlist[[c("TFoverrepresented")]]<- TFtable$overRepresented
  }
  returnerlist[[c("summary")]]<-summary
  return(returnerlist)
  
}


###plot_cluster_profiler information
#cluster_name<-c("turquoise")
plot_single_cluster<-function(igraph_object=igraph_object,
                              cluster_name=c("all"),                 #,"cluster_name"
                              top_percentage_for_hubs=0.25,
                              allowed_edges=25,
                              allowed_edges_between_hubs=1,
                              string_needs_to_be_redone=TRUE,
                              string_treshold=400,
                              color_STRING=c("#c95555"),
                              color_edges=c("grey"),
                              no_strings_to_be_string_hub=1,
                              label_all_TF=TRUE,
                              color_label_if_TF=c("#D3436EFF"),
                              color_label_normal=c("#231151FF"),
                              width_label_edge=2,
                              width_normal=1,
                              percentage_named=0.1,
                              size_label_boxes=20){
  summary<-allargs()
  if(length(cluster_name)>1){
    cluster_name<-paste0(cluster_name,collapse = "|")
    cluster_genes<-cluster_information[grepl(cluster_name,cluster_information$try),c("Gene")]
  }else{
    cluster_genes<-cluster_information[grepl(cluster_name,cluster_information$try),c("Gene")]
  }
  
  edges<-correlation_df[correlation_df$correlation>chosen_cutoff,]
  edges<-edges[edges$from %in% cluster_genes ,]
  edges<-edges[edges$to %in% cluster_genes ,]
  
  edges<-edges[base::order(edges$correlation),]
  edges$rank<-1:nrow(edges)
  
  igraph_object<-igraph::graph_from_data_frame(edges,vertices = NULL,directed=FALSE)
  #cluster network before
  # e <- igraph::get.edgelist(igraph_object,names = FALSE)
  # l <- qgraph.layout.fruchtermanreingold(e,vcount=igraph::vcount(igraph_object),weights = igraph::edge.attributes(igraph_object)$weight,
  #                                        area=15*(igraph::vcount(igraph_object)^2*2),repulse.rad=(igraph::vcount(igraph_object)^3.1),
  #                                        niter = 100)
  # plot(igraph_object,layout=l,vertex.label=igraph::V(igraph_object)$name,vertex.frame.color="black",vertex.size=5,
  #      vertex.color= ifelse(igraph::V(igraph_object)$name %in% TF_list$Mouse,"red","floralwhite"))
  # dev.off()
  # 
  
  #layout=layout_nicely(igraph_object)
  degrees<-as.data.frame(igraph::degree(igraph_object,mode=c("all"),normalized = FALSE))
  degrees$Gene<-rownames(degrees)
  degrees$cumulataive_correlation<-apply(degrees,1,function(x){
    cumulataive_correlation_df<-edges[(edges$from %in% x[2] |edges$to %in% x[2]) ,] 
    sum(cumulataive_correlation_df$correlation)/nrow(cumulataive_correlation_df)
  })
  colnames(degrees)[1]<-c("degrees")
  hub_genes_1<-degrees[degrees$degrees>(1-top_percentage_for_hubs)*max(degrees$degrees) & degrees$cumulataive_correlation>(1-top_percentage_for_hubs)*max(degrees$cumulataive_correlation),]
  #hub_genes_2<-degrees[degrees$degrees>(1-(top_percentage_for_hubs+0.15))*max(degrees$degrees) & degrees$cumulataive_correlation>(1-(top_percentage_for_hubs+0.15))*max(degrees$cumulataive_correlation),]
  #hub_genes_2<-hub_genes_2[!(hub_genes_2$Gene %in% (hub_genes_1$Gene)),]
  
  # rest<-degrees$Gene[!(degrees$Gene%in%hub_genes_1$Gene |degrees$Gene%in%hub_genes_2$Gene )]
  rest<-degrees$Gene[!(degrees$Gene%in%hub_genes_1$Gene )]
  rest<-degrees[rest,]
  
  edges_selected<-list()
  
  for(hub_genes_1_genes in hub_genes_1$Gene){
    
    hub_genes_1_edges<-edges[(grepl(hub_genes_1_genes,edges$from)|grepl(hub_genes_1_genes,edges$to)),]
    hub_genes_1_edges<-hub_genes_1_edges[base::order(hub_genes_1_edges$correlation,decreasing = TRUE),]
    tmp<-hub_genes_1_edges[1:allowed_edges,]
    genes_involved<-as.data.frame(unique(c(as.character(tmp$from),as.character(tmp$to))))
    genes_involved<-as.data.frame(genes_involved[!is.na(genes_involved)])
    rownames(genes_involved)<-genes_involved[,1]
    genes_involved<-genes_involved[!(genes_involved[,1]%in%hub_genes_1_genes),]
    genes_duble<-as.character(genes_involved[genes_involved %in% hub_genes_1$Gene])
    edges_simple<-hub_genes_1_edges[1:allowed_edges,][!(hub_genes_1_edges[1:allowed_edges,]$from%in%genes_duble |hub_genes_1_edges[1:allowed_edges,]$to%in%genes_duble  ),]
    edges_simple_add<-hub_genes_1_edges[1:allowed_edges_between_hubs,][(hub_genes_1_edges[1:allowed_edges_between_hubs,]$from%in%genes_duble |hub_genes_1_edges[1:allowed_edges_between_hubs,]$to%in%genes_duble  ),]
    edges_simple$from<-as.character(edges_simple$from)
    edges_simple$to<-as.character(edges_simple$to)
    edges_selected[[hub_genes_1_genes]]<-edges_simple
    edges_simple_add$from<-as.character(edges_simple_add$from)
    edges_simple_add$to<-as.character(edges_simple_add$to)
    edges_selected[[paste0(hub_genes_1_genes,"_add")]]<-edges_simple_add
  }
  
  # for(hub_genes_1_genes in hub_genes_2$Gene){
  #   hub_genes_1_edges<-edges[(grepl(hub_genes_1_genes,edges$from)|grepl(hub_genes_1_genes,edges$to)),]
  #   hub_genes_1_edges<-hub_genes_1_edges[base::order(hub_genes_1_edges$correlation,decreasing = TRUE),]
  #   tmp<-hub_genes_1_edges[1:10,]
  #   genes_involved<-as.data.frame(unique(c(as.character(tmp$from),as.character(tmp$to))))
  #   genes_involved<-as.data.frame(genes_involved[!is.na(genes_involved)])
  #   rownames(genes_involved)<-genes_involved[,1]
  #   genes_involved<-genes_involved[!(genes_involved[,1]%in%hub_genes_1_genes),]
  #   genes_duble<-as.character(genes_involved[genes_involved %in% hub_genes_2$Gene])
  #   edges_simple<-hub_genes_1_edges[1:10,][!(hub_genes_1_edges[1:10,]$from%in%genes_duble|hub_genes_1_edges[1:10,]$to%in%genes_duble ),]
  #   edges_simple_add<-hub_genes_1_edges[1,][(hub_genes_1_edges[1,]$from%in%hub_genes_1$Gene |hub_genes_1_edges[1,]$to%in%hub_genes_1$Gene  ),]
  #   edges_simple$from<-as.character(edges_simple$from)
  #   edges_simple$to<-as.character(edges_simple$to)
  #   edges_selected[[hub_genes_1_genes]]<-edges_simple
  #   edges_simple_add$from<-as.character(edges_simple_add$from)
  #   edges_simple_add$to<-as.character(edges_simple_add$to)
  #   edges_selected[[paste0(hub_genes_1_genes,"_add")]]<-edges_simple_add
  # }
  
  for(hub_genes_1_genes in rest$Gene){
    hub_genes_1_edges<-edges[(grepl(hub_genes_1_genes,edges$from)|grepl(hub_genes_1_genes,edges$to)),]
    hub_genes_1_edges<-hub_genes_1_edges[base::order(hub_genes_1_edges$correlation,decreasing = TRUE),]
    tmp<-hub_genes_1_edges[1:15,]
    genes_involved<-as.data.frame(unique(c(as.character(tmp$from),as.character(tmp$to))))
    genes_involved<-as.data.frame(genes_involved[!is.na(genes_involved)])
    
    rownames(genes_involved)<-genes_involved[,1]
    genes_involved<-genes_involved[!(genes_involved[,1]%in%hub_genes_1_genes),]
    genes_duble<-as.character(genes_involved[genes_involved %in% rest$Gene])
    edges_simple<-hub_genes_1_edges[1,][!(hub_genes_1_edges[1,]$from%in%genes_duble|hub_genes_1_edges[1,]$to%in%genes_duble ),]
    # if(nrow(edges_simple)==0){
    #   edges_simple$
    #   edges_selected[[hub_genes_1_genes]]<-edges_simple
    # }else{
    edges_selected[[hub_genes_1_genes]]<-edges_simple
    #}
  }
  
  edges_selected<-list.rbind(edges_selected)
  
  edges_selected<-unique.data.frame(edges_selected)
  rownames(edges_selected)<-1:nrow(edges_selected)
  
  edges_selected$index<-1:nrow(edges_selected)
  unique_genes<-unique(c(edges_selected$from , edges_selected$to))
  unique_genes<-unique_genes[unique_genes%in%rest$Gene]
  list_new_edges<-list()
  
  for(every_gene in unique_genes){
    tmp_1<-edges_selected[grepl(every_gene,edges_selected$from),]
    tmp_2<-edges_selected[grepl(every_gene,edges_selected$to),]
    tmp<-rbind(tmp_1,tmp_2)
    tmp<-tmp[order(tmp$correlation,decreasing = TRUE),]
    list_new_edges[[every_gene]]<-tmp[-(1),]
  }
  
  delted_edges<-list.rbind(list_new_edges)
  edges_selected<-edges_selected[-(delted_edges$index),]
  edges_selected$index<-NULL
  edges_selected$col<-color_edges
  ###############################################################STRING INPUT##############################################
  if(string_needs_to_be_redone ==TRUE){
    edges_selected_rest<-list()
    for(hub_genes_1_genes in rest$Gene){
      hub_genes_1_edges<-edges[(grepl(hub_genes_1_genes,edges$from)|grepl(hub_genes_1_genes,edges$to)),]
      hub_genes_1_edges<-hub_genes_1_edges[base::order(hub_genes_1_edges$correlation,decreasing = TRUE),]
      
      edges_selected_rest[[hub_genes_1_genes]]<-hub_genes_1_edges
      
    }
    
    # for(hub_genes_1_genes in hub_genes_1$Gene){
    #   hub_genes_1_edges<-edges[(grepl(hub_genes_1_genes,edges$from)|grepl(hub_genes_1_genes,edges$to)),]
    #   hub_genes_1_edges<-hub_genes_1_edges[base::order(hub_genes_1_edges$correlation,decreasing = TRUE),]
    #   
    #   edges_selected_rest[[hub_genes_1_genes]]<-hub_genes_1_edges
    #   
    # }
    edges_rest<- unique(list.rbind(edges_selected_rest))
    ##
    #nodesvec <- as.character(unique(unique(edges_rest$from) , unique(edges_rest$to)))
    nodesvec<-V(igraph_object)$name
    library(STRINGdb)
    if(organism=="Mouse"){
      string_db <- STRINGdb$new(version = "10" , species = 10090, score_threshold=0.7)
    }else{
      string_db <- STRINGdb$new(version = "10" , species = 9606, score_threshold=0.7)
    }
    node <- as.data.frame(nodesvec)
    colnames(node) <- "Gene"
    example1_mapped <- string_db$map( node, "Gene", removeUnmappedRows = TRUE )
    interactions <- string_db$get_interactions(example1_mapped$STRING_id)
    ####
    #interactions_test <- string_db$get_link(example1_mapped$STRING_id)
    ####
    interactions$from<-apply(interactions,1,function(x){
      if(x[1] %in% example1_mapped$STRING_id){
        example1_mapped[grepl(paste0("^",x[1],"$"),example1_mapped$STRING_id),][1,1]
      }
    })
    interactions$to<-apply(interactions,1,function(x){
      if(x[2] %in% example1_mapped$STRING_id){
        example1_mapped[grepl(paste0("^",x[2],"$"),example1_mapped$STRING_id),][1,1]
      }
    })
    #######set threshold score to be an link 
    interactions<-as.data.frame(interactions[interactions$combined_score >string_treshold,])
    
    ######
    clust_fr_to <- data.frame(get.edgelist(igraph_object),stringsAsFactors = FALSE)
    colnames(clust_fr_to) <- c("from", "to")
    clust_fr_to$col <- c("")
    clust_fr_to$string <- paste(clust_fr_to$from, clust_fr_to$to , sep = " ")
    #clust_fr_to$string_rev <- paste(clust_fr_to$to, clust_fr_to$from , sep = " ")
    interactions$string <- paste(interactions$from, interactions$to , sep = " ")
    interactions$string_rev<-paste(interactions$to, interactions$from , sep = " ")
    ####
    
    clust_fr_to$col<-apply(clust_fr_to,1,function(x){
      if(x[4] %in% interactions$string | x[4] %in% interactions$string_rev){
        color_STRING
      }else{
        color_edges
      }
    })
    
   # FromTo<-merge(edges_rest,clust_fr_to, by=c("from","to"))
    FromTo<-clust_fr_to
    #####################ACVHTUNG
    #FALL ABKLÄREN WENN KEINE KNOWNS
    if(color_STRING %in% FromTo$col){
      FromTo<-FromTo[grepl(color_STRING,FromTo$col),]
      
      FromTo$string<-NULL
      
      
      ######find out how often single gene appears
      gene_conn_known<-as.data.frame(table(c(as.character(FromTo$from), as.character(FromTo$to))))
      gene_conn_known<-gene_conn_known[order(gene_conn_known$Freq,decreasing=TRUE),]
      
      keep_weil_hub_genes<-gene_conn_known[(gene_conn_known$Freq> no_strings_to_be_string_hub),] ### änern zu absolute zahl!
      keep_weil_hub<-FromTo[(FromTo$from %in% keep_weil_hub_genes$Var1 | FromTo$to %in% keep_weil_hub_genes$Var1),]
      ##delete weil hub string to hub string
      #keep_weil_hub<-keep_weil_hub[!(keep_weil_hub$from %in% keep_weil_hub_genes$Var1 & keep_weil_hub$to %in% keep_weil_hub_genes$Var1),]
      
      
      keep_weil_con_to_known<-FromTo[(FromTo$from %in% hub_genes_1$Gene | FromTo$to %in% hub_genes_1$Gene ),]#|FromTo$from %in% hub_genes_2$Gene | FromTo$to %in% hub_genes_2$Gene ),]
      
      
      keep_from_string<-rbind(keep_weil_hub,keep_weil_con_to_known)
      keep_from_string<-unique.data.frame(keep_from_string)
      keep_from_string<-data.frame(keep_from_string[,1:2],correlation=1,
                                   pval=0,
                                   rank=1,
                                   col=keep_from_string$col)
      
      
      edges_selected<-rbind(edges_selected,keep_from_string)
      
    }else{
      print(paste0("String has found no known interactions above certain threshold: ",string_treshold))
    } 
  }
  
  ###################################################################################################################
  edges_selected<- edges_selected[complete.cases(edges_selected),]
  
  edges_selected$index<-1:nrow(edges_selected)
  edges_selected_doubles<-edges_selected[duplicated.data.frame(edges_selected[,1:2]),]
  
  to_delete<-list()
  if(nrow(edges_selected_doubles)>0){
    for(double_edges in 1:nrow(edges_selected_doubles)){
      to_check<-c(edges_selected_doubles[double_edges,c("from")],edges_selected_doubles[double_edges,c("to")])
      edges_doubeld<-edges_selected[edges_selected$from %in%to_check & edges_selected$to %in%to_check,]
      tmp<-edges_doubeld[grep(color_edges,edges_doubeld$col),]
      edges_selected<-edges_selected[-(grep(tmp$index,edges_selected$index)),]
    }
  }
  
  
  ########################
  
  igraph_object<-igraph::graph_from_data_frame(edges_selected,vertices = NULL,directed=FALSE)
  degree_table<-data.frame(degree=igraph::degree(igraph_object),gene=igraph::V(igraph_object)$name)
  
  edges_selected$weight<-0.5
  # edges_selected$weight<-apply(edges_selected,1,function(x){
  #  if(x[1] %in% rest$Gene | x[2] %in% rest$Gene){
  #    100
  #  }else{
  #    if((x[1] %in% hub_genes_1$Gene &  x[2] %in% hub_genes_2$Gene) | (x[1] %in% hub_genes_2$Gene & x[2] %in% hub_genes_1$Gene)){
  #      80
  #    }else{
  #      0.001
  #    }
  #  }
  #   }
  # )
  
  edges_selected$weight<-apply(edges_selected,1,function(x){
    if(x[1] %in% rest$Gene | x[2] %in% rest$Gene){
      1
    }else{
      
      0.001
    }
  }
  )
  
  edges_selected$curved<-apply(edges_selected,1,function(x){
    if(x["col"]==color_edges){
      TRUE  
    }else{
      FALSE
    }
  })
  
  igraph_object<-igraph::graph_from_data_frame(edges_selected,vertices = NULL,directed=FALSE)
  
  
  
  #layout<-layout_with_fr(igraph_object ,niter = 10000)
  vertex_attributes<-data.frame(Gene=igraph::V(igraph_object)$name)
  igraph_object<-igraph::delete.vertices(igraph_object,is.na(igraph::V(igraph_object)$name))
  
  degrees<-as.data.frame(igraph::degree(igraph_object,mode=c("all"),normalized = FALSE))
  
  degrees$Gene<-igraph::V(igraph_object)$name
  
  degrees$bins <- as.numeric(cut(degrees[,1], breaks=fields::stats.bin(1:max(degrees[,1]),degrees[,1],N=10)$breaks,include.lowest = TRUE))+0.5
  
  colnames(degrees)<-c("degree","Gene","bins")
  
  
  vertex_attributes<-merge(vertex_attributes,degrees,by="Gene")
  rownames(vertex_attributes)<-vertex_attributes$Gene
  vertex_attributes<-vertex_attributes[c(igraph::V(igraph_object)$name),]
  #kleine zahl = kleinste degree
  # vertex_attributes$vertex_size<-apply(vertex_attributes,1,function(x){
  #   if(x[1] %in% hub_genes_1$Gene){
  #     10
  #   }else{
  #     if(x[1] %in% hub_genes_2$Gene){
  #       5
  #     }else{
  #       2
  #     }
  #   }
  # }
  # )
  
  
  vertex_attributes$color<-"white"
  
  
  ####label df
  vertex_attributes$label=c("FALSE")
  vertex_attributes$shape<-"circle"
  nodes_added<-list()
  edges_added<-list()
  no<-0
  
  if(label_all_TF==TRUE){
    for(Gene in as.character(vertex_attributes$Gene)){
      tmp<-vertex_attributes[grepl(paste0("^",Gene,"$"),vertex_attributes$Gene),]
      if(tmp$Gene %in% TF_list[,organism]){ # vorher 8
        nodes_added[[tmp$Gene]]<-data.frame(row.names=paste0(tmp$Gene,"_label_weil_TF"),degree=1,Gene=tmp$Gene,bins=max(vertex_attributes$bins),color=c("white"),label=c("TRUE"),shape="rectangle")
        no<-no+1
        edges_added[[tmp$Gene]]<-data.frame(from=paste0(tmp$Gene,"_label_weil_TF"),to=tmp$Gene, correlation=1, pval=0.05,rank=0,col=c("black"),index=max(edges_selected$index)+no,weight=2,curved=c("FALSE")) #group=c("weight")
      }
    }
  }
  
  for(Gene in as.character(vertex_attributes$Gene)){
    tmp<-vertex_attributes[grepl(paste0("^",Gene,"$"),vertex_attributes$Gene),]
    if(tmp$bins > max(vertex_attributes$bins)*(1-percentage_named)){ # vorher 8
      nodes_added[[tmp$Gene]]<-data.frame(row.names=paste0(tmp$Gene,"_label"),degree=1,Gene=tmp$Gene,bins=max(vertex_attributes$bins),color=c("white"),label=c("TRUE"),shape="rectangle")
      no<-no+1
      edges_added[[tmp$Gene]]<-data.frame(from=paste0(tmp$Gene,"_label"),to=tmp$Gene, correlation=1, pval=0.05,rank=0,col=c("black"),index=max(edges_selected$index)+no,weight=2,curved=c("FALSE")) #,group=c("weight")
    }
  }
  
  
  
  nodes_added<-list.rbind(nodes_added)
  
  
  edges_added<-list.rbind(edges_added)
  
  edges_added$col<-apply(edges_added,1,function(x){
    if(x["to"] %in% TF_list[,organism]){
      color_label_if_TF
    }else{
      color_label_normal
    }
  })
  
  edges_added$edge_width<-width_label_edge
  
  edges_selected$edge_width<-width_normal
  
  
  vertex_attributes_label<-rbind(vertex_attributes,nodes_added)
  edges_selected_label<-rbind(edges_selected,edges_added)
  
  
  vertex_attributes_label$label<-apply(vertex_attributes_label,1,function(x){
    if(x["label"]==c("TRUE")){
      x[1]
    }else{
      c("")
    }
    
  })
  
  plot(igraph_object)
  
  vertex_attributes_label$size<-0
  vertex_attributes_label$size2<-0
  vertex_attributes_label[vertex_attributes_label$shape == "rectangle",]$size<-apply(vertex_attributes_label[vertex_attributes_label$shape == "rectangle",],1,function(x){
    (strwidth (x[1]) + strwidth ("oo") )*size_label_boxes
  })
  vertex_attributes_label[vertex_attributes_label$shape == "rectangle",]$size2<-strwidth ("M") * 2*size_label_boxes
  vertex_attributes_label[vertex_attributes_label$shape =="circle",]$size<-as.numeric(vertex_attributes_label[vertex_attributes_label$shape =="circle",]$bins) 
  vertex_attributes_label[vertex_attributes_label$shape =="circle",]$size2<-as.numeric(vertex_attributes_label[vertex_attributes_label$shape =="circle",]$bins) 
  
  
  igraph_object<-igraph::graph_from_data_frame(edges_selected_label,vertices = NULL,directed=FALSE)
  
  
  vertex_attributes_label<-vertex_attributes_label[igraph::V(igraph_object)$name,]
  
  
  
  vertex_attributes_label$label_color<-apply(vertex_attributes_label,1,function(x){
    if(x["label"] %in% TF_list[,organism]){
      color_label_if_TF
    }else{
      color_label_normal
    }
  })
  
  vertex_attributes_label$label_TF<-apply(vertex_attributes_label,1,function(x){
    if(x["Gene"] %in% as.character(vertex_attributes_label[duplicated(vertex_attributes_label$Gene),"Gene"] ) &   x["shape"]==c("circle")){
      c("TRUE")
    }else{
      c("FALSE")
    }
  })
  
  
  vertex_attributes_label$frame_color<-apply(vertex_attributes_label,1,function(x){
    if((x["label"] %in% TF_list[,organism] & x["shape"]==c("rectangle")) | (x["Gene"] %in% TF_list[,organism] & x["label_TF"]==c("TRUE") )){
      color_label_if_TF
    }else{
      if(x["shape"]==c("rectangle") | x["label_TF"]==c("TRUE") ){
        color_label_normal
      }else{
        c("grey26")}
    }
  })
  
  HUBS<-unique(cluster_genes)
  
  Hubs<-HUBS
  #HUBS<-as.character(vertex_attributes_label[vertex_attributes_label$shape==c("rectangle"),c("Gene")])
  ##########################HIERACHYY
  Hub_Hoods_df<-Hierarchy_df(Hubs,chosen_cutoff)
  Hub_Hoods_df<-Hub_Hoods_df$BayesDataFrame
  
  databox<-boxplot(Hub_Hoods_df$Ranking)
  
  Hub_Hoods_df[Hub_Hoods_df$Ranking>databox$stats[5,1],"Ranking"]<-databox$stats[5,1]
  Hub_Hoods_df[Hub_Hoods_df$Ranking<databox$stats[1,1],"Ranking"]<-databox$stats[1,1]
  
  ###
  color_hierachy<-ggplot2::ggplot(Hub_Hoods_df,aes(x=Hub_Hoods_df$No.parents,y=Hub_Hoods_df$No.children))+geom_point(aes(color=Hub_Hoods_df$Ranking))+
    scale_colour_gradient2(low=c("#5A85E1"),mid=c("white"),high=c("#792187"),
                          guide="colourbar")
  
  gg <- ggplot_build(color_hierachy)

  farben_from_hierachy<-unlist(gg$data[[1]]["colour"],use.names = FALSE)
  

  # rbPal <- colorRampPalette(c('#5A85E1','#92AADD','#B1B1B2','#C873D6','#792187'))
  # Hub_Hoods_df$colour<-pal(10)[as.numeric(cut(Hub_Hoods_df$Ranking,breaks = 10))]
  Hub_Hoods_df$colour<-farben_from_hierachy
  
  Genes_to_color<-unique(as.character(vertex_attributes_label[!(vertex_attributes_label$frame_color=="grey26"),"Gene"]))
  Hub_Hoods_df_to_use<-Hub_Hoods_df[as.character(Hub_Hoods_df$Gene) %in% Genes_to_color,]
  
  #########binning data into 10 bins
  # minimum_ranking<-min(Hub_Hoods_df_to_use$Ranking)
  # maximum_ranking<-max(Hub_Hoods_df_to_use$Ranking)
  # 
  # bin_it(minimum_ranking,maximum_ranking,10)
  
  #rbPal <- colorRampPalette(c('#5A85E1','#92AADD','#B1B1B2','#C873D6','#792187'), bias =3)
  
  
  #Hub_Hoods_df2$colour <- rbPal(50000)[as.numeric(cut(Hub_Hoods_df2$Ranking,breaks = 50000))]
  ###################################################################################################################
  Hub_Hoods_df_to_use<-as.data.frame(cbind(as.character(Hub_Hoods_df_to_use$Gene),as.character(Hub_Hoods_df_to_use$colour)))
  
  
  ##ACHTUNG NEU
  # vertex_attributes_label$color_from_hierachy<-apply(vertex_attributes_label,1,function(x){
  #   if(x["label_TF"]==c("TRUE")){
  #     tmp<-Hub_Hoods_df_to_use[grepl(x["Gene"],Hub_Hoods_df_to_use$V1),]
  #     as.character(tmp$V2)
  #   }else{
  #     if(x["shape"]==c("rectangle")){
  #       c("#F8F8FF") #weiß
  #     }else{
  #     c("#B1B1B2") #grau
  #     }
  #   }
  # })
  
  vertex_attributes_label$color_from_hierachy<-apply(vertex_attributes_label,1,function(x){
    if(x["shape"]==c("circle") & x["Gene"] %in%Hub_Hoods_df_to_use$V1 ){
      tmp<-Hub_Hoods_df_to_use[grepl(paste0("^",as.character(x["Gene"]),"$"),as.character(Hub_Hoods_df_to_use$V1)),]
      as.character(tmp$V2)
    }else{
      if(x["shape"]==c("rectangle")){
        c("#F8F8FF")
      }else{
        c("#B1B1B2")
      }
    }
  })
  
  
  
  
  vertex_attributes_label<-vertex_attributes_label[igraph::V(igraph_object),]
  list_of_results<-list()
  list_of_results[[c("igraph_object_small")]]<-igraph_object
  list_of_results[["vertex_attributes_cluster"]]<-vertex_attributes_label
  list_of_results[["edges_attributes"]]<-edges_selected_label
  list_of_results[["summary"]]<-summary
  
  return(list_of_results)
  
  # library(qgraph)
  # e <- get.edgelist(igraph_object,names = FALSE)
  # l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(igraph_object),weights = edge.attributes(igraph_object)$weight,
  #                                        area=10*(vcount(igraph_object)^2*2),repulse.rad=(vcount(igraph_object)^3.1),
  #                                        niter = 10000)
  # 
  # 
  # 
  # 
  # # plot(0, type="n", ann=FALSE, axes=FALSE, xlim=extendrange(l[,1]), 
  # #      ylim=extendrange(l[,2]))
  # plot(igraph_object,layout=l,
  #      vertex.size=vertex_attributes_label$size,
  #      vertex.size2=vertex_attributes_label$size2,
  #      vertex.label=vertex_attributes_label$label,
  #      vertex.color=vertex_attributes_label$color_from_hierachy,
  #      vertex.label.color=vertex_attributes_label$label_color,
  #      vertex.frame.color=vertex_attributes_label$frame_color,
  #      edge.color=edges_selected_label$col,edge.curved=edges_selected_label$curved,
  #      vertex.label.font=2,vertex.shape=vertex_attributes_label$shape,
  #      edge.width=edges_selected_label$edge_width)
  # 
  dev.off()
}







############plot cluster profiler information

cluster_information_to_color<-function(vertex_attributes=list_of_results[[c("vertex_attributes_cluster")]],
                                       clusterprofiler_results=clusterprofiler_results,
                                       colors_for_anno=data.frame(type=c(names(clusterprofiler_results)),color=c("green",
                                                                                                        "orange",
                                                                                                        "blue",
                                                                                                        "yellow",
                                                                                                        "purple",
                                                                                                        "red"))){
  summary<-allargs()
  color_saved<-list()
  for(data_frame_info in names(clusterprofiler_results)){
    if(!(data_frame_info==c("TFoverrepresented"))){
    data_frame_tmp<-list()
    for(rows in 1:nrow(clusterprofiler_results[[data_frame_info]])){
      Genes<-strsplit(clusterprofiler_results[[data_frame_info]]$geneID[rows],"/")[[1]]
      name<-clusterprofiler_results[[data_frame_info]]$Description[rows]
      data_frame_tmp[[name]]<-data.frame(Genes,type=rep(name,length.out=length(Genes)))
    }
    
    data_frame_tmp<-list.rbind(data_frame_tmp)
    
    
    color_saved[[data_frame_info]]<-apply(vertex_attributes,1,function(x){
      if(x["Gene"] %in% as.character(data_frame_tmp$Genes)){
        res<-data_frame_tmp[grepl(paste0("^",x["Gene"],"$"),data_frame_tmp$Genes),]
        if(nrow(res)>1){
          paste0(res$type, collapse = "/")
        }else{
          as.character(res$type)
        }
      }else{
        c(" ")
      }
    }
    )
    
    
    tmp<-data.frame(color_saved[[data_frame_info]],stringsAsFactors = FALSE)
    
    color_saved[[paste0(data_frame_info,"_Color")]]<-apply(tmp,1,function(x){
      if(!(x[1])==c(" ")){
        as.character(colors_for_anno[grepl(data_frame_info,colors_for_anno$type),2])
      }else{
        c("grey")
      }
    }
    )
  }
  }
  
  color_saved<-list.cbind(color_saved)
  rownames(color_saved)<-vertex_attributes$Gene
  RT<-list()
  RT[["summary"]]<-summary
  RT[["color_saved"]]<-color_saved
  return(RT)
}



# Hallmark_df<-list()
# for(rows in 1:nrow(clusterprofiler_results$Hallmark)){
#   Genes<-strsplit(clusterprofiler_results$Hallmark$geneID[rows],"/")[[1]]
#   name<-clusterprofiler_results$Hallmark$ID[rows]
#   Hallmark_df[[name]]<-data.frame(Genes,type=rep(name,length.out=length(Genes)))
# }
# 
# Hallmark_df<-list.rbind(Hallmark_df)
# 
# 
# vertex_attributes$color_hallmark<-apply(vertex_attributes,1,function(x){
#   if(x[1] %in% Hallmark_df$Genes){
#     res<-Hallmark_df[grepl(paste0("^",x[1],"$"),Hallmark_df$Genes),]
#     if(nrow(res)>1){
#       paste0(res$type, collapse = "/")
#     }else{
#       as.character(res$type)
#     }
#     }else{
#       c("white")
#   }
# }
# )


##Maria PPi
# library(STRINGdb)
# EdgeColPPI <- function(FromTo) {
#   
#   
#   
#   nodesvec <- unique(unique(FromTo$from) , unique(FromTo$to))
#   #library(STRINGdb)
#   if(organism=="Mouse"){
#     string_db <- STRINGdb$new(version = "10" , species = 10090)
#   }else{
#     string_db <- STRINGdb$new(version = "10" , species = 9606)
#   }
#  
#   node <- as.data.frame(nodesvec)
#   colnames(node) <- "Gene"
#   example1_mapped <- string_db$map( node, "Gene", removeUnmappedRows = TRUE )
#   interactions <- string_db$get_interactions(example1_mapped$STRING_id)
#   interactions$from<-apply(interactions,1,function(x){
#     if(x[1] %in% example1_mapped$STRING_id){
#       example1_mapped[grepl(x[1],example1_mapped$STRING_id),][1,1]
#     }
#   })
#   interactions$to<-apply(interactions,1,function(x){
#     if(x[2] %in% example1_mapped$STRING_id){
#       example1_mapped[grepl(x[2],example1_mapped$STRING_id),][1,1]
#     }
#   })
#   clust_fr_to <- cbind.data.frame(FromTo$from , FromTo$to)
#   colnames(clust_fr_to) <- c("from", "to")
#   clust_fr_to$col <- c("")
#   clust_fr_to$string <- paste(clust_fr_to$from, clust_fr_to$to , sep = " ")
#   interactions$string <- paste(interactions$from, interactions$to , sep = " ")
#   clust_fr_to$col<-apply(clust_fr_to,1,function(x){
#     if(x[4] %in% interactions$string){
#       c("orange")
#     }else{
#       c("grey")
#     }
#   })
#   
#   FromTo<-merge(FromTo,clust_fr_to, by=c("from","to"))
#   
#   FromTo <- FromTo[,!names(FromTo)=="string"]
#   
#   
#   return(FromTo)
#   
#   
# }



#####base_area
#HUBS === VEKTOR !!!

Hierarchy <- function(Hubs , chosen_cutoff=chosen_cutoff){
  library(bnlearn)
  cutted_data <- correlation_df[(correlation_df$correlation>=chosen_cutoff),]
  Hub_1_df <- cutted_data[(cutted_data$from %in% Hubs | cutted_data$to %in% Hubs ),]
  Hub_1_df <- Hub_1_df[complete.cases(Hub_1_df),]
  Hub_connections <- as.vector(unique(unique(Hub_1_df$from), unique(Hub_1_df$to)))
  Hub_expression <- as.data.frame(original_data[(rownames(original_data) %in% Hub_connections),])
  Hub_expression <- as.data.frame(t(Hub_expression))
  bnnet<-hc(Hub_expression)
  netlist <- bnnet$nodes
  Hub_Hoods <- netlist[list.names(netlist) %in% Hubs]
  Hub_Hoods_df <- as.data.frame(list.names(Hub_Hoods))
  colnames(Hub_Hoods_df) <- "Gene"
  flat <- list.flatten(Hub_Hoods)
  a <- 1:length(flat)
  b <- a[seq(3, length(a) , 4)]
  no <- 0
  liste <- list()
  for(x in b) {
    no <- no+1
    liste[[no]]<- paste(length(flat[[x]]))}
  No.parents <- as.vector(as.numeric(liste[1:length(Hubs)]))
  Hub_Hoods_df$No.parents <- No.parents
  a <- 1:length(flat)
  b <- a[seq(4, length(a) , 4)]
  no <- 0
  liste <- list()
  for(x in b) {
    no <- no+1
    liste[[no]]<- paste(length(flat[[x]]))}
  No.children <- as.vector(as.numeric(liste[1:length(Hubs)]))
  Hub_Hoods_df$No.children <- No.children
  Hub_Hoods_df$No.neighbours <- Hub_Hoods_df$No.children+Hub_Hoods_df$No.parents
  Hub_Hoods_df$Ranking <- (Hub_Hoods_df$No.children/(Hub_Hoods_df$No.parents+1))*(Hub_Hoods_df$No.neighbours/max(Hub_Hoods_df$No.neighbours)) #um zu vermeiden, dass der Nenner null ist
  rbPal <- colorRampPalette(c('#B1B1B2','#982D80FF'))
  Hub_Hoods_df$colour <- rbPal(50)[as.numeric(cut(Hub_Hoods_df$Ranking,breaks = 50))]
  return(as.data.frame(Hub_Hoods_df))
}

# library(viridis)
# rev(magma(8))[c(-1,-8)]

 

summary_wrapped_up<-function(print_text=TRUE){
  
  tmp<-data.frame(original_no_genes= nrow(Dataset_1),
                                 original_no_samples = ncol(Dataset_1),
                                 correlation_measure = summary$correlation$type_of_correlation,
                                 correlation_cutoff = chosen_cutoff,
                                 no_of_nodes_after_cutoff=length(levels(summary$plot_network$data$from)),
                                 no_of_edges_after_cutoff=nrow(summary$plot_network$data),
                                 plotting_layout_chosen=summary$plot_network$layout,
                                 condition_used_for_GFC=summary$GFC_calculation$group,
                                 setted_range_for_GFC=summary$GFC_calculation$range_GFC,
                                 network_cluster_algorithm=summary$heatmap_clustered$cluster_algo,
                                 specific_cluster_chosen=summary$`I-GIN`$cluster_name,
                                 IGIN_no_of_nodes=length(V(summary$Output_IGIN)),
                                 IGIN_no_of_edges=length(E(summary$Output_IGIN)))
  
  
  if(print_text==TRUE){
    print(paste0("Original Data had ", tmp$original_no_genes, " Genes and " , tmp$original_no_samples, " Samples"))
    print(paste0(tmp$correlation_measure, " Correlation measure was used, found correlations were kept if their correlation was over ",tmp$correlation_cutoff))
    print(paste0("After cutting, ", tmp$no_of_nodes_after_cutoff, " Genes were left and ", tmp$no_of_edges_after_cutoff , " Edges"))
    print(paste0("Network was plotted with ", tmp$plotting_layout_chosen ," layout and clustered with ", tmp$network_cluster_algorithm))
    print(paste0("The condition ",tmp$condition_used_for_GFC, " is basis of the GFC-calculations" ))
    print(paste0("Cluster ", tmp$specific_cluster_chosen, " was chosen to be investigated further"))
    print(paste0("Resulting I-GIN plot existed out of ", tmp$IGIN_no_of_nodes, " Genes and ",tmp$IGIN_no_of_edges, " Edges"))
  }   
  summary_df<-t(tmp)
  colnames(summary_df)<-c("setted_values")
  
  return(summary_df)
}





#####search and label certain genes

add_gene_label<-function(gene_to_add=c("")){
  vertex_attributes_change=c("FALSE")
  for(gene in gene_to_add){
  if(gene %in% vertex_attributes_label$Gene ){
    tmp<-vertex_attributes_label[grepl(paste0("^",gene,"$"),vertex_attributes_label$Gene),]
    tmp$Gene<-as.character(tmp$Gene)
    if(nrow(tmp)>1){
      print(paste0(gene," is already labeled - have another look at the plot"))#
      vertex_attributes_change=c("TRUE")
    }else{
      vertex_attributes_label[vertex_attributes_label$Gene==gene,c("label")]<-gene
      vertex_attributes_label[vertex_attributes_label$Gene==gene,c("label_color")]<-"red"
      vertex_attributes_change=c("TRUE")
    }
  }else{
    print(paste0(gene," is not in the plot at all, if it was in the cluster and you want it to be included, try loosen parameters in function:plot_single_cluster"))
  }
  }
  if(vertex_attributes_change==c("TRUE")){
    plot(igraph_object_small,layout=l,
         vertex.size=vertex_attributes_label$size,
         vertex.size2=vertex_attributes_label$size2,
         vertex.label=vertex_attributes_label$label,
         vertex.color=vertex_attributes_label$color_from_hierachy,
         vertex.label.color=vertex_attributes_label$label_color,
         vertex.frame.color=vertex_attributes_label$frame_color,
         edge.color=edges_selected_label$col,edge.curved=edges_selected_label$curved,
         vertex.label.font=2,vertex.shape=vertex_attributes_label$shape,
         edge.width=edges_selected_label$edge_width,
         main=to_save["cluster_name",])
  }
}




#############################################################################################################################################
#############################################################################################################################################
###############################################   MARIE WORK  ###############################################################################
#############################################################################################################################################
#############################################################################################################################################
#Hubs<-as.character(unique(cluster_information[cluster_information$try=="turquoise","Gene"]))

Hierarchy_df <- function(Hubs, chosen_cutoff){
  library(bnlearn)
  library(rlist)
  cutted_data <- correlation_df[(correlation_df$correlation>=chosen_cutoff),]
  
  Hub_1_df2 <- cutted_data[(cutted_data$from %in% Hubs & cutted_data$to %in% Hubs ),]
  Hub_1_df2 <- Hub_1_df2[complete.cases(Hub_1_df2),]
  Hub_connections2 <- data.frame(Gene=unique(c(unique(as.character(Hub_1_df2$from)),unique(as.character(Hub_1_df2$to)))))
  Hub_connections2 <- as.vector(Hub_connections2$Gene)
  Hub_expression2 <- as.data.frame(original_data[(rownames(original_data) %in% Hub_connections2),])
  Hub_expression2 <- as.data.frame(t(Hub_expression2))
  #bnnet2<-hc(Hub_expression2)
  #bnnet2 <- mmhc(Hub_expression2)
  bnnet2 <- bnlearn::hc(Hub_expression2)
  # class(bnnet)
  
  # shor
  #bnnet2 <- bnlearn::mmpc(Hub_expression2)
  netlist2 <- bnnet2$nodes
  Hub_Hoods2 <- netlist2[list.names(netlist2) %in% Hubs]
  Hub_Hoods_df2 <- as.data.frame(list.names(Hub_Hoods2))
  colnames(Hub_Hoods_df2) <- "Gene"
  flat2 <- list.flatten(Hub_Hoods2)
  a <- 1:length(flat2)
  b <- a[seq(3, length(a) , 4)]
  no <- 0
  liste2 <- list()
  for(x in b) {
    no <- no+1
    liste2[[no]]<- paste(length(flat2[[x]]))}
  No.parents2 <- as.vector(as.numeric(liste2[1:length(liste2)]))
  Hub_Hoods_df2$No.parents <- No.parents2
  a <- 1:length(flat2)
  b <- a[seq(4, length(a) , 4)]
  no <- 0
  liste2 <- list()
  for(x in b) {
    no <- no+1
    liste2[[no]]<- paste(length(flat2[[x]]))}
  No.children2 <- as.vector(as.numeric(liste2[1:length(liste2)]))
  Hub_Hoods_df2$No.children <- No.children2
  Hub_Hoods_df2$No.neighbours <- Hub_Hoods_df2$No.children+Hub_Hoods_df2$No.parents
  
  Hub_Hoods_df2$Ranking <- (Hub_Hoods_df2$No.children-(Hub_Hoods_df2$No.parents))*(Hub_Hoods_df2$No.children/(Hub_Hoods_df2$No.parents+1)) #um zu vermeiden, dass der Nenner null ist
  
  #Hub_Hoods_df2[Hub_Hoods_df2$Ranking<0,"Ranking"]<-(Hub_Hoods_df2$No.children-(Hub_Hoods_df2$No.parents))*(Hub_Hoods_df2$No.parents/(Hub_Hoods_df2$No.children+1))

  Hub_Hoods_df2$Ranking<-apply(Hub_Hoods_df2,1,function(x){
    if(as.numeric(x["Ranking"])<0){
      (as.numeric(x["No.children"])-as.numeric(x["No.parents"]))*(as.numeric(x["No.parents"])/(as.numeric(x["No.children"])+1))
    }else{
      as.numeric(x["Ranking"])
    }
  })
  
  ##scaling ranking
  
  
  rbPal <- colorRampPalette(c('#5A85E1','#92AADD','#B1B1B2','#C873D6','#792187'), bias =3)
  Hub_Hoods_df2$colour <- rbPal(100)[as.numeric(cut(Hub_Hoods_df2$Ranking,breaks = 100))]
  #return(as.data.frame(Hub_Hoods_df2))
  BayesListAndDataframe <- list()
  BayesListAndDataframe[["BayesDataFrame"]]<- as.data.frame(Hub_Hoods_df2)
  BayesListAndDataframe[["BayesCompleteList"]]<- flat2
  return(BayesListAndDataframe)
  #return(bnnet2)
}



arcs <- function(BayesListe) {
  c <- 1:length(BayesListe)
  d <- c[seq(4, length(c) , 4)]
  
  lieblibngs_liste<-list()
  no<-0
  
  
  for(listenplatz in d){
    
    no<-no+1
    tmp<-BayesListe[[listenplatz]]
    if(length(tmp)==0){
      print(listenplatz)
    }else{
      name<-names( BayesListe[listenplatz])
      dataframe_yey<-data.frame(from=rep(name,length(tmp )),to=as.character(tmp))
      lieblibngs_liste[[no]]<-dataframe_yey
    }
  }
  Marie_FROM_TO<-list.rbind(lieblibngs_liste)  
  
  Marie_FROM_TO$from <- gsub(pattern = ".children" , replacement = "" , Marie_FROM_TO$from)
  
  from_to <- Marie_FROM_TO
  
  return(from_to)}






arcs2 <- function(BayesListe) {
  c <- 1:length(BayesListe)
  d <- c[seq(3, length(c) , 4)]
  
  lieblibngs_liste<-list()
  no<-0
  
  
  for(listenplatz in d){
    
    no<-no+1
    tmp<-BayesListe[[listenplatz]]
    if(length(tmp)==0){
      print(listenplatz)
    }else{
      name<-names( BayesListe[listenplatz])
      dataframe_yey<-data.frame(from=rep(name,length(tmp )),to=as.character(tmp))
      lieblibngs_liste[[no]]<-dataframe_yey
    }
  }
  Marie_FROM_TO<-list.rbind(lieblibngs_liste)  
  
  Marie_FROM_TO$from <- gsub(pattern = ".parents" , replacement = "" , Marie_FROM_TO$from)
  
  from_to <- Marie_FROM_TO
  
  return(from_to)}





circosPlot <- function(clust.No, cutoff , CPratio , PCratio){
  
  library(circlize)
  
  Hubs <-as.vector(cluster_information[cluster_information$`clp$membership`==clust.No, 4])
  
  #Hub_Hoods_df2<-Hierarchy_df(Hubs = Hubs, chosen_cutoff = cutoff)
  
  #flat2 <- Hierarchy_list(Hubs = Hubs, chosen_cutoff = cutoff)
  
  BayesListAndDataframe <-Hierarchy_df(Hubs = Hubs, chosen_cutoff = cutoff)
  
  Hub_Hoods_df2 <- BayesListAndDataframe[[1]]
  
  flat2 <- BayesListAndDataframe[[2]]
  
  from_to <- arcs(flat2) 
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(Hub_Hoods_df2$Ranking,decreasing=T)),]
  
  Hub_Hoods_df2$Group <- "0"
  
  vec1 <- as.vector(Hub_Hoods_df2[(Hub_Hoods_df2$No.children/Hub_Hoods_df2$No.parents)>=CPratio,1])
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(Hub_Hoods_df2$No.parents,decreasing=F)),]
  
  vecP <- as.vector(Hub_Hoods_df2[Hub_Hoods_df2$No.parents==0,1])
  
  Hub_Hoods_df2$Group <- apply(Hub_Hoods_df2,1, function(x){
    if(x[1] %in% vec1){
      Hub_Hoods_df2$Group <- as.character(match(x[1], vec1))
    }else{
      "0"}})
  
  
  
  
  
  
  liste <- list()
  no <- 0
  for(a in vec1){
    no <- no+1
    ft <- from_to[from_to$from==a ,]
    tmp <- Hub_Hoods_df2[Hub_Hoods_df2$Gene %in% ft$to,]
    tmp$from <- a
    tmp$Group <- as.character(tmp$Group)
    tmp$Group <- apply(tmp , 1 , function(x){
      if(x[7]== "0"){
        as.character(match(as.character(x[8]), vec1))
        #"T"
      }else{
        x["Group"]
        #"F"
      }})
    
    liste[[no]] <- tmp
  }
  
  yay<-as.data.frame(list.rbind(liste))
  
  lost_genes <- Hub_Hoods_df2[Hub_Hoods_df2$Gene %in%vec1,]
  
  lost_genes$from <- "x"
  
  lost_genes$Gene <- as.character(lost_genes$Gene)
  
  yay <-rbind(yay, lost_genes)
  
  #dup <- yay[duplicated(yay$Gene)==TRUE,]
  
  #not_dup <-yay[duplicated(yay$Gene)==FALSE,]
  
  
  
  
  nochnelist <- list()
  no <-0
  for( a in unique(as.character(yay$Gene))){
    no <- no+1
    tmp <- yay[yay$Gene==a,]
    tmp <- tmp[tmp$Group == min(tmp$Group),]
    nochnelist[[no]] <- tmp
    
  }
  
  yay2 <-as.data.frame(list.rbind(nochnelist))
  
  yay2 <- yay2[,1:7]
  
  final_df <- unique(yay2)
  
  Hub_Hoods_df2 <- final_df
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(as.numeric(Hub_Hoods_df2$Group),decreasing=T)),]
  
  farbvector <- if(length(vec1)==1){c("#70014b")}else{if(length(vec1)<=10){
    as.vector(rainbow(length(vec1)))
  }else{
    colorRampPalette(c("yellow","red", "green","blue","orange","pink","black"))(length(vec1))
  }}
  
  farb_df<-data.frame(Gene=vec1, Colour=farbvector)
  
  Hub_Hoods_df2$colour <- apply(Hub_Hoods_df2,1,function(x){
    if(x[1] %in% farb_df$Gene){
      Hub_Hoods_df2$colour <- as.character(farb_df[farb_df$Gene==x[1],2])
    }else{
      Hub_Hoods_df2$colour <- "lightgrey"}})
  
  H_H_df2<-Hub_Hoods_df2[complete.cases(Hub_Hoods_df2),]
  
  new_from_to <- from_to[(from_to$from %in% H_H_df2$Gene & from_to$to %in% H_H_df2$Gene),]
  
  n_n_f_t <- new_from_to[(new_from_to$from %in% vec1 ),]
  
  vertex <- as.vector(unique(n_n_f_t$from))
  vertex2<-as.vector(unique(n_n_f_t$to))
  vertex <- unique(c(vertex,vertex2))
  H <- H_H_df2[(H_H_df2$Gene %in% vertex),]
  
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE )
  g <- chordDiagram(n_n_f_t, order = H$Gene , grid.col = H$colour , transparency = 0.55 , self.link = 1, symmetric = T, annotationTrack = "grid", preAllocateTracks = 1, link.lty = if(length(vec1)==1){1}else{0} , link.border = "#660043")
  g <- circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, cex=0.5 , facing = "clockwise", niceFacing = TRUE, adj = c(-0.6,-0.1))
    circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA )
  g <- title(paste("Cluster ",as.character(clust.No),": ","Controlling Hot Spots of Bayesian Network and their Children")) 
  
  #pdf(g,"circosPlotParents.pdf", paper = "a4")
  
  
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[,1:6]
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(Hub_Hoods_df2$No.parents,decreasing=F)),]
  
  vecP <- as.vector(Hub_Hoods_df2[Hub_Hoods_df2$No.parents==0,1])
  
  
  Hub_Hoods_df2$Group <- apply(Hub_Hoods_df2,1, function(x){
    if(x[1] %in% vecP){
      Hub_Hoods_df2$Group <- as.character(match(x[1], vecP))
    }else{
      "0"}})
  
  
  
  
  
  
  liste <- list()
  no <- 0
  for(a in vecP){
    no <- no+1
    ft <- from_to[from_to$from==a ,]
    tmp <- Hub_Hoods_df2[Hub_Hoods_df2$Gene %in% ft$to,]
    tmp$from <- a
    tmp$Group <- as.character(tmp$Group)
    tmp$Group <- apply(tmp , 1 , function(x){
      if(x[7]== "0"){
        as.character(match(as.character(x[8]), vecP))
        #"T"
      }else{
        x["Group"]
        #"F"
      }})
    
    liste[[no]] <- tmp
  }
  
  yay<-as.data.frame(list.rbind(liste))
  
  lost_genes <- Hub_Hoods_df2[Hub_Hoods_df2$Gene %in%vecP,]
  
  lost_genes$from <- "x"
  
  lost_genes$Gene <- as.character(lost_genes$Gene)
  
  yay <-rbind(yay, lost_genes)
  
  #dup <- yay[duplicated(yay$Gene)==TRUE,]
  
  #not_dup <-yay[duplicated(yay$Gene)==FALSE,]
  
  
  
  
  nochnelist <- list()
  no <-0
  for( a in unique(as.character(yay$Gene))){
    no <- no+1
    tmp <- yay[yay$Gene==a,]
    tmp <- tmp[tmp$Group == min(tmp$Group),]
    nochnelist[[no]] <- tmp
    
  }
  
  yay2 <-as.data.frame(list.rbind(nochnelist))
  
  yay2 <- yay2[,1:7]
  
  final_df <- unique(yay2)
  
  Hub_Hoods_df2 <- final_df
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(as.numeric(Hub_Hoods_df2$Group),decreasing=T)),]
  
  farbvector <- if(length(vecP)==1){c("#70014b")}else{if(length(vecP)<=10){
    as.vector(rainbow(length(vecP)))
  }else{
    colorRampPalette(c("yellow","red", "green","blue","pink","orange","black"))(length(vecP))
  }}
  
  farb_df<-data.frame(Gene=vecP, Colour=farbvector)
  
  Hub_Hoods_df2$colour <- apply(Hub_Hoods_df2,1,function(x){
    if(x[1] %in% farb_df$Gene){
      Hub_Hoods_df2$colour <- as.character(farb_df[farb_df$Gene==x[1],2])
    }else{
      Hub_Hoods_df2$colour <- "lightgrey"}})
  
  H_H_df2<-Hub_Hoods_df2[complete.cases(Hub_Hoods_df2),]
  
  new_from_to <- from_to[(from_to$from %in% H_H_df2$Gene & from_to$to %in% H_H_df2$Gene),]
  
  n_n_f_t <- new_from_to[(new_from_to$from %in% vecP ),]
  
  vertex <- as.vector(unique(n_n_f_t$from))
  vertex2<-as.vector(unique(n_n_f_t$to))
  vertex <- unique(c(vertex,vertex2))
  H <- H_H_df2[(H_H_df2$Gene %in% vertex),]
  
  circos.clear()
  circos.par(start.degree = 90, clock.wise = TRUE )
  h <- chordDiagram(n_n_f_t, order = H$Gene , grid.col = H$colour , transparency = 0.55 , self.link = 1, symmetric = T, annotationTrack = "grid", preAllocateTracks = 1, link.lty = if(length(vecP)==1){1}else{0} , link.border = "#660043")
  h <- circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, cex=0.7 , facing = "clockwise", niceFacing = TRUE, adj = c(-0.6,-0.1))
    circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA )
  h <- title(paste("Cluster ",as.character(clust.No),": ","Roots of Bayesian Network and their Children")) 
  
  #pdf(g,"circosPlotParents.pdf", paper = "a4")
  
  from_to <- arcs2(flat2) 
  
  Hub_Hoods_df2 <- Hub_Hoods_df2[(order(Hub_Hoods_df2$Ranking,decreasing=F)),]
  
  Hub_Hoods_df2$Group <- "0"
  
  vecReceiver <- as.vector(Hub_Hoods_df2[(Hub_Hoods_df2$No.parents/Hub_Hoods_df2$No.children)>=PCratio,1])
  
  BayesListAndDataframe[["GenesHotSpots"]] <-vec1
  
  BayesListAndDataframe[["GenesParents"]] <- vecP
  
  BayesListAndDataframe[["GenesReceiver"]] <- vecReceiver
  
  return(BayesListAndDataframe)
  
  
  
}


netPlot <- function(data , CPratio , PCratio){
  
  library(rlist)
  library(bnlearn)
  
  c<-as.vector(as.numeric(data[as.numeric(length(data)),]))
  
  net <- bnlearn::hc(data)
  
  Hub_Hoods2 <- net$nodes
  
  Hub_Hoods_df2 <- as.data.frame(list.names(Hub_Hoods2))
  
  colnames(Hub_Hoods_df2) <- "Gene"
  
  flat2 <- list.flatten(Hub_Hoods2)
  
  a <- 1:length(flat2)
  b <- a[seq(3, length(a) , 4)]
  no <- 0
  liste2 <- list()
  for(x in b) {
    no <- no+1
    liste2[[no]]<- paste(length(flat2[[x]]))}
  No.parents2 <- as.vector(as.numeric(liste2[1:length(liste2)]))
  Hub_Hoods_df2$No.parents <- No.parents2
  a <- 1:length(flat2)
  b <- a[seq(4, length(a) , 4)]
  no <- 0
  liste2 <- list()
  for(x in b) {
    no <- no+1
    liste2[[no]]<- paste(length(flat2[[x]]))}
  
  No.children2 <- as.vector(as.numeric(liste2[1:length(liste2)]))
  Hub_Hoods_df2$No.children <- No.children2
  Hub_Hoods_df2$No.neighbours <- Hub_Hoods_df2$No.children+Hub_Hoods_df2$No.parents
  Hub_Hoods_df2$Ranking <- (Hub_Hoods_df2$No.children/(Hub_Hoods_df2$No.parents+1))*(Hub_Hoods_df2$No.neighbours/max(Hub_Hoods_df2$No.neighbours))
  rbPal <- colorRampPalette(c('#FFFFFF','#982D80FF'))
  Hub_Hoods_df2 <- Hub_Hoods_df2[order(Hub_Hoods_df2$Ranking, decreasing=T),]
  rownames(Hub_Hoods_df2) <- 1:nrow(Hub_Hoods_df2)
  Hub_Hoods_df2$colour <- rbPal(50)[as.numeric(cut(Hub_Hoods_df2$Ranking,breaks = 50))]
  BayesListAndDataframe <- list()
  
  BayesListAndDataframe[["BayesDataFrame"]]<- as.data.frame(Hub_Hoods_df2)
  
  BayesListAndDataframe[["BayesCompleteList"]] <- net
  
  vec1 <- as.vector(Hub_Hoods_df2[(Hub_Hoods_df2$No.children/Hub_Hoods_df2$No.parents)>=CPratio,1])
  
  vecP <- as.vector(Hub_Hoods_df2[Hub_Hoods_df2$No.parents==0,1])
  
  vecReceiver <- as.vector(Hub_Hoods_df2[(Hub_Hoods_df2$No.parents/Hub_Hoods_df2$No.children)>=PCratio,1])
  
  BayesListAndDataframe[["GenesHotSpots"]] <-vec1
  
  BayesListAndDataframe[["GenesParents"]] <- vecP
  
  BayesListAndDataframe[["GenesReceiver"]] <- vecReceiver
  
  igraph<-graph_from_data_frame(d=net$arcs , vertices = Hub_Hoods_df2$Gene)
  
  plot.new()
  
  Hub_Hoods_df2$size <- apply(Hub_Hoods_df2,1, function(x){
    (strwidth(x[1]) + strwidth("oo",font=2) )* 2*50})
  
  plot.igraph(igraph, vertex.color = Hub_Hoods_df2$colour, 
              vertex.label.color = "#2d0024", 
              vertex.shape = "rectangle", vertex.size = Hub_Hoods_df2$size,
              edge.arrow.width = 1 , edge.arrow.size =0.4, vertex.label.family="Arial", edge.width=2,
              edge.curved=0)
  
  BayesOutput <- BayesListAndDataframe
  
  return(BayesOutput)
}



BayesNet <- function(clust.name , cutoff , CPratio , PCratio) {
  
  library(bnstruct)
  
  
  cfg_df <- as.data.frame(cluster_information[cluster_information$try == clust.name,4])
  
  colnames(cfg_df) <- "Gene"
  
  nodes <- as.vector(rownames(original_data[rownames(original_data) %in% cfg_df$Gene, ]))
  
  test <- as.data.frame(original_data[rownames(original_data) %in% cfg_df$Gene, ])
  
  test2 <- as.data.frame(t(test))
  
  test2["max",] <- apply(test2[,],2,max)
  
  no.nodes <- length(nodes)
  
  if (no.nodes < 70) {
    
    BayesOutput <- netPlot(data = test2 , CPratio = CPratio , PCratio = PCratio)
    
    return(BayesOutput)
    
  }else{
    
    clust.no<-unique(cluster_information[grepl(paste0("^",clust.name,"$"),cluster_information$try),1])
    
    BayesOutput <- circosPlot(clust.No = clust.no, cutoff = cutoff, CPratio, PCratio)
    
    return(BayesOutput)
    
  }
  
}


#BayesOutput <- BayesNet(clust.no = 6 , cutoff = chosen_cutoff , CPratio = 5 , PCratio = 3) 


###########################
bayesian_based_subnetworks<-function(igraph_object=igraph_object_small){
  selected_gene_parent<-as.character(vertex_attributes_label$Gene[vertex_attributes_label$color_from_hierachy=="#982D80"])
  print(paste0(length(selected_gene_parent)," Nodes are selceted as top influencer based on bayesian approach"))
  list_of_subnetworks<-make_ego_graph(igraph_object,order=1, nodes = selected_gene_parent)
  
  ##corresponding GFC's
  rbPal <- colorRampPalette(c('#0000FF','#FFFFFF','#FFFFFF','#FF0000'))
  tmp_GFC<-GFC_all_genes
  tmp_GFC$Gene<-NULL
  
  corresponding_GFC_Color<-list()
  for(list_element in 1:length(list_of_subnetworks)){
    tmp_list<-list()
    tmp_network<-list_of_subnetworks[[list_element]]
    tmp_genes<-V(tmp_network)$name
    tmp_genes<-tmp_genes[!(grepl("label",tmp_genes))]
    for(tmp_gene in tmp_genes){
      GFC_value<-GFC_all_genes[tmp_gene,]
      GFC_value$Gene<-NULL
      GFC_color<-rbPal(8)[cut(as.numeric(GFC_value),breaks=c(min(tmp_GFC),-2,-1.5,-1,0,1,2,1.5,max(tmp_GFC)),include.lowest = TRUE)]
      tmp_list[[tmp_gene]]<-GFC_color
    }
    corresponding_GFC_Color[[list_element]]<-tmp_list
  }
  values_for_pie_piece_size<-list()
  for(no_of_subnetworks in 1:length(corresponding_GFC_Color)){
    values_for_pie_piece_size[[no_of_subnetworks]]<-rep((1/ncol(tmp_GFC)*10),times=ncol(tmp_GFC))
    
    
  }
  
  
  
  par(mfrow=c(2,2))
  
  for(i in 1:length(list_of_subnetworks))  {
    to_be_delted_nodes<-V(list_of_subnetworks[[i]])$name[grepl("label",V(list_of_subnetworks[[i]])$name)]
    if(length(to_be_delted_nodes)>0){
      tmp_edgelist<-as.data.frame(get.edgelist(list_of_subnetworks[[i]]))
      tmp_edgelist<-tmp_edgelist[!(grepl(to_be_delted_nodes,tmp_edgelist$V1)),]
      tmp_edgelist<-tmp_edgelist[!(grepl(to_be_delted_nodes,tmp_edgelist$V2)),]
      tmp_edgelist<-as.matrix(tmp_edgelist)
      tmp_network<-graph_from_edgelist(tmp_edgelist,directed = FALSE)
    }else{
      tmp_network<-list_of_subnetworks[[i]]
    }
    
    
    
    # if(length(to_be_delted_nodes)>0){
    #   for(each in to_be_delted_nodes){
    #     to_delete_1<-paste0(each,"|",strsplit(each,"_")[[1]][1])
    #     to_delete_2<-paste0(strsplit(each,"_")[[1]][1],"|",each)
    #     igraph::delete.edges(list_of_subnetworks[[i]],to_delete_1)
    #     igraph::delete.edges(list_of_subnetworks[[i]],to_delete_2)
    #     tmp_network<-igraph::delete.vertices(list_of_subnetworks[[i]],each)
    #   }
    # }
    
    E(tmp_network)$curved=FALSE
    for(vertex in V(tmp_network)$name){
      V(tmp_network)[vertex]$pie.color<-corresponding_GFC_Color[[i]][vertex]
    }
    
    plot(tmp_network,
         layout=layout_as_tree,
         vertex.label.color="white",
         vertex.label=V(tmp_network)$name,
         vertex.shape="pie",
         vertex.pie=values_for_pie_piece_size,
         vertex.size=50,
         main=as.character(i))
    
  }
}

####
merge_cluster<-function(data=clustered_heatmap_data,cluster_wanted=2){
  summary<-allargs()
  new <- kmeans(data, centers = cluster_wanted)
  new_cluster <- new["cluster"]
  new_cluster <- as.data.frame(new_cluster[["cluster"]])
  colnames(new_cluster) <- "new_cluster"
  new_cluster$old <- rownames(new_cluster)
  
  data$cluster <- rownames(data) 
  data <- data[order(data$cluster),]
  new_cluster <- new_cluster[order(new_cluster$old),]  
  data$new <- new_cluster$new_cluster
  data <- data[order(data$new),]
  
  list1 <- list()
  
  for(n in 1:length(unique(data$new))){
    list1[[n]] <- data[data$new==n,1:4]
  }
  
  
  list2 <- list()
  
  
  for(a in 1:length(list1)){
    
    tmp <- list1[[a]]
    newtmp <- apply(tmp,2,mean)
    newtmp <- t(newtmp)
    list2[[a]]<- newtmp
    
  }
  
  newmat <- as.data.frame(list.rbind(list2))
  
  rnames <- data[,c("cluster","new")]
  rnames$x <- str_extract(rnames$cluster, "\\[\\d+\\]")
  rnames$size <- as.numeric(str_extract_all(rnames$x, "\\d+"))
  
  no<-0
  liste0<-list()
  for(x in 1:as.numeric(max(rnames$new))){
    no<-no+1
    
    liste0[[no]]<-sum(rnames[rnames$new==x,4])}
  
  rnames$x<- apply(rnames,1,function(a){ 
    if(as.numeric(a[2]) %in% 1:as.numeric(max(rnames$new))){
      a[3] <- liste0[[as.numeric(a[2])]]
    }})
  
  rnames$merged <- paste("cluster ", rnames$new,  " : ", rnames$x, " genes")
  
  
  names <- unique(rnames$merged)
  
  rownames(newmat) <- names
  
  p<-pheatmap::pheatmap(newmat, cluster_rows = T, cluster_cols = T, color= col.pal, border_color = "grey")
  
  rnames$old <- sapply(strsplit(rownames(rnames),"\\["),"[",1)
  rnames$cl.colour<-"c"
  col_vec <- labels2colors(unique(rnames$new), zeroIsGrey = T)
  rnames$cl.colour <- apply(rnames,1,function(a){ 
    if(as.numeric(a[2]) %in% 1:as.numeric(max(rnames$new))){
      as.character(col_vec[as.numeric(a[2])])
    }})
  
  cluster_information$try2 <- "c"
  
  cluster_information$try2 <- apply(cluster_information,1,function(a){
    
    if(a[2] %in% rnames$old){
      as.character(rnames[rnames$old==a[2],7])
    }else{a[5]}
  }
  )
  
  cluster_information$clp2<-"c"
  
  cluster_information$clp2 <- apply(cluster_information,1,function(a){
    
    if(a[5] %in% rnames$cl.colour){
      as.character(unique(rnames[rnames$cl.colour==a[5],2]))
    }else{a[6]}
  }
  )
  
  cluster_information$`clp$membership`<- cluster_information$clp2 # Anmerkung: c könnten im weiteren Verlauf ein Problem werden -> Zeilen löschen?
  
  cluster_information$try <- cluster_information$try2 
  
  
  rnames$merged <- paste("cluster ", rnames$new, " ",rnames$cl.colour, " : ", rnames$x, " genes")
  
  
  names <- unique(rnames$merged)
  
  rownames(newmat) <- names
  
  p<-pheatmap::pheatmap(newmat, cluster_rows = T, cluster_cols = T, color= col.pal, border_color = "grey")
  
  cluster_information$try2<-NULL
  cluster_information$clp2<-NULL
  
  cluster_information[cluster_information$`clp$membership`==c("c"),c("try")]<-c("grey")
  cluster_information[cluster_information$`clp$membership`==c("c"),c("clp$membership")]<-0
  
  RT<-list()
  RT[["cluster_information_new"]]<-cluster_information
  RT[["summary"]]<-summary
  return(RT)
}



require(plyr)
merge_cluster_Thomas<-function(data=clustered_heatmap_data, height = height, cluster_cols = T){

  clusters <- hclust(dist(clustered_heatmap_data))
  

  new_cluster <- dendextend:::cutree.dendrogram(as.dendrogram(clusters),h=height)
  
  
  new_cluster <- as.data.frame(new_cluster)
  colnames(new_cluster) <- "new_cluster"
  new_cluster$old <- rownames(new_cluster)
  
  data$cluster <- rownames(data) 
  data <- data[order(data$cluster),]
  new_cluster <- new_cluster[order(new_cluster$old),]  
  data$new <- new_cluster$new_cluster
  data <- data[order(data$new),]
  
  list1 <- list()
  
  for(n in 1:length(unique(data$new))){
    list1[[n]] <- data[data$new==n,1:(ncol(data)-2)]
  }
  
  
  list2 <- list()
  
  
  for(a in 1:length(list1)){
    
    tmp <- list1[[a]]
    newtmp <- apply(tmp,2,mean)
    newtmp <- t(newtmp)
    list2[[a]]<- newtmp
    
  }
  
  newmat <- as.data.frame(list.rbind(list2))
  
  rnames <- data[,c("cluster","new")]
  rnames$x <- str_extract(rnames$cluster, "\\[\\d+\\]")
  rnames$size <- as.numeric(str_extract_all(rnames$x, "\\d+"))
  
  no<-0
  liste0<-list()
  for(x in 1:as.numeric(max(rnames$new))){
    no<-no+1
    
    liste0[[no]]<-sum(rnames[rnames$new==x,4])}
  
  rnames$x<- apply(rnames,1,function(a){ 
    if(as.numeric(a[2]) %in% 1:as.numeric(max(rnames$new))){
      a[3] <- liste0[[as.numeric(a[2])]]
    }})
  
  rnames$merged <- paste("cluster ", rnames$new,  " : ", rnames$x, " genes")
  
  
  names <- unique(rnames$merged)
  
  rownames(newmat) <- names
  
  #p<-pheatmap::pheatmap(newmat, cluster_rows = T, cluster_cols = T, color= col.pal, border_color = "grey")
  
  rnames$old <- sapply(strsplit(rownames(rnames),"\\["),"[",1)
  rnames$cl.colour<-"c"
  col_vec <- labels2colors(unique(rnames$new), zeroIsGrey = T)
  rnames$cl.colour <- apply(rnames,1,function(a){ 
    if(as.numeric(a[2]) %in% 1:as.numeric(max(rnames$new))){
      as.character(col_vec[as.numeric(a[2])])
    }})
  
  
  #cluster_information$try2 <- "c"
  
  
  cluster_information$`clp$membership` <- apply(cluster_information,1,function(a){
    
    if(a[2] %in% rnames$old){
      as.character(rnames[rnames$old==a[2],2])
    }else{a[5]}
  }
  )
  
  if("waste" %in% cluster_information$try){
  cluster_information[cluster_information$try=="waste",]$`clp$membership` <- 0
  }
  
  cluster_information$try2 <- apply(cluster_information,1,function(a){
    
    if(a[1] %in% rnames$new){
      as.character(rnames[rnames$old==a[2],7])
    }else{a[5]}
  }
  )
  
  if("white" %in% cluster_information$vertex_color){
  cluster_information[cluster_information$vertex_color=="white",]$try2 <- "white"
  }
  # cluster_information$`clp$membership` <- apply(cluster_information,1,function(a){
  #   
  #   if(a[2] %in% rnames$old){
  #     as.character(rnames[rnames$old==a[2],2])
  #   }else{a[5]}
  # }
  # )
  
  cluster_information$clp2<-"c"
  
  cluster_information$clp2 <- apply(cluster_information,1,function(a){
    
    if(a[5] %in% rnames$cl.colour){
      as.character(unique(rnames[rnames$cl.colour==a[5],2]))
    }else{a[6]}
  }
  )
  
  #cluster_information$`clp$membership`<- cluster_information$clp2 # Anmerkung: c könnten im weiteren Verlauf ein Problem werden -> Zeilen löschen?
  
  cluster_information$try <- cluster_information$try2 
  cluster_information$vertex_color_old <- cluster_information$vertex_color 
  cluster_information$vertex_color <- cluster_information$try2 
  
  
  rnames$merged <- paste("cluster ", rnames$new, " ",rnames$cl.colour, " : ", rnames$x, " genes")
  
  
  names <- unique(rnames$merged)
  
  rownames(newmat) <- names
  
  p<-pheatmap::pheatmap(newmat, cluster_rows = T, cluster_cols = cluster_cols, color= col.pal, border_color = "grey")
  
  
  cluster_information$try2<-NULL
  cluster_information$clp2<-NULL
  
  
 
  
  par(mfrow=c(1,2))
  frame()

  heatmap<-Heatmap(newmat, 
                   cluster_rows = T,
                   cluster_columns = T,
                   color = col.pal)
  gb = grid.grabExpr(draw(heatmap),newpage = FALSE)#color auch hier ersetzt
  #upViewport()
  grid.arrange(gb,ncol=2)
  
  Sys.sleep(5)
  plot(igraph_object,layout=layout,
       vertex.label=NA,
       edge.color="lightgrey",
       vertex.shape="fcircle",
       vertex.size=2,
       vertex.color=cluster_information$vertex_color,
       vertex.frame.color="lightgrey")
  
  
  
  RT<-list()
  RT[["cluster_information_new"]]<-cluster_information
  RT[["summary"]]<-summary
  RT[["plot"]]<-p
  
  return(RT)
  
}




compareClusterGO <- function(){
  
  cl_info_compare <- cluster_information[!cluster_information$try=="waste",]
  entrez_all = clusterProfiler::bitr(cl_info_compare$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  entrez_all$CLUSTER <- apply(entrez_all,1,function(x){
    
    as.character(cl_info_compare[cl_info_compare$Gene==x[1],"try"])
    
  })
  entrez_list <- list()
  for(x in entrez_all$CLUSTER){
    
    tmp <- as.vector(entrez_all[entrez_all$CLUSTER==x,"ENTREZID"])
    tmp <- as.character(tmp) 
    entrez_list[[x]] <- tmp
    
  }
  
  
  library(clusterProfiler)
  #data(gcSample)
  y=compareCluster(entrez_list, fun='enrichGO', pvalueCutoff=0.01, pAdjustMethod="BH", OrgDb=org.Hs.eg.db,ont="BP",readable=T)
  dotplot(y, showCategory=3)
  
  y_compare <- y@compareClusterResult
  
  return(y_compare)
}





clusterCircos <- function(hm_vec, seg, range, y_comp){
  
  GO_immune_response <- unique(GO_immune_response$Annotated.Term)
  GO_metabolic_process <- unique(GO_metabolic_process$Annotated.Term)
  GO_signalling <- unique(GO_signalling$Annotated.Term)
  
  number_TF <- list()
  number_Epi <- list()
  
  cluster_information <- cluster_information[!cluster_information$try=="white",]
  
  for(x in unique(cluster_information$try)){
    
    tmp <- cluster_information[cluster_information$try == x,"Gene"]
    #tmp <- tmp$Gene
    TFs <- tmp[tmp %in%  TF_list$Human]
    Epis <- tmp[tmp %in% epigenetic_modulators$gene]
    number_TF[[x]] <- c(length(TFs), length(tmp))
    number_Epi[[x]] <- c(length(Epis), length(tmp))
    
  }
  
  number_TF <- as.data.frame(list.rbind(number_TF))
  colnames(number_TF) <- c("part","total")
  number_TF$x <- rownames(number_TF)
  number_TF$names <- paste0(number_TF$x," ","[",number_TF$total,"]")
  number_Epi <- as.data.frame(list.rbind(number_Epi))
  colnames(number_Epi) <- c("part","total")
  number_Epi$x <- rownames(number_Epi)
  number_Epi$names <- paste0(number_Epi$x," ","[",number_Epi$total,"]")
  
  f2 = factor(unique(cluster_information$try))
  
  GFC_list <- list()
  
  for(x in unique(cluster_information$try)){
    
    tmp <- cluster_information[cluster_information$try == x,"Gene"]
    #tmp <- tmp$Gene
    GFC <- as.data.frame(GFC_all_genes[rownames(GFC_all_genes) %in% tmp,1:(ncol(GFC_all_genes)-1)])
    GFC <- apply(GFC,2,mean)
    GFC_list[[x]] <- GFC
  }
  
  GFC <- as.data.frame(list.rbind(GFC_list))
  
  if(hm_vec == TRUE){
    
    GFC <- GFC[heatmap_vec]
    
  }else{
    
    
  }
  
  
  df <- data.frame(x=f2, y=number_TF, z = number_Epi)
  
  df <- df[sort(levels(f2)),]
  
  GFC <- as.data.frame(t(GFC))
  
  GFC <- GFC[levels(df$x)]
  
  colnames(GFC) <- df$y.names
  
  GFC <- as.data.frame(t(GFC))
  
  df$g <- GFC
  
  vec <- c(seq(2,24,3))
  
  bed2 <- list()
  
  for(a in rownames(GFC)){
    
    tmp <- number_Epi[number_Epi$names==a,1:2]
    m <- as.data.frame(matrix(c(rep(1,tmp[1,"part"]), rep(0,max(number_Epi$part)-tmp[1,"part"])), ncol = 1))
    n <- as.data.frame(matrix(c(rep(a,tmp[1,"part"]), rep(a,max(number_Epi$part)-tmp[1,"part"])), ncol = 1))
    p <- as.data.frame(matrix(c(rep(0,tmp[1,"part"]), rep(0,max(number_Epi$part)-tmp[1,"part"])), ncol = 1))
    q <- as.data.frame(matrix(c(rep(1,tmp[1,"part"]), rep(1,max(number_Epi$part)-tmp[1,"part"])), ncol = 1))
    m <- cbind(n,p,q,m)
    bed2[[a]] <- m
    
  }
  
  
  bed2 <- list.rbind(bed2)
  rownames(bed2) <- 1:nrow(bed2)
  colnames(bed2) <- c("chr", "start", "end", "value1")
  
  bed3 <- list()
  
  for(a in rownames(GFC)){
    
    tmp <- number_TF[number_TF$names==a,1:2]
    m <- as.data.frame(matrix(c(rep(1,tmp[1,"part"]), rep(0,max(number_TF$part)-tmp[1,"part"])), ncol = 1))
    n <- as.data.frame(matrix(c(rep(a,tmp[1,"part"]), rep(a,max(number_TF$part)-tmp[1,"part"])), ncol = 1))
    p <- as.data.frame(matrix(c(rep(0,tmp[1,"part"]), rep(0,max(number_TF$part)-tmp[1,"part"])), ncol = 1))
    q <- as.data.frame(matrix(c(rep(1,tmp[1,"part"]), rep(1,max(number_TF$part)-tmp[1,"part"])), ncol = 1))
    m <- cbind(n,p,q,m)
    bed3[[a]] <- m
    
  }
  
  
  bed3 <- list.rbind(bed3)
  rownames(bed3) <- 1:nrow(bed3)
  colnames(bed3) <- c("chr", "start", "end", "value1")
  
  
  
  bed4 <- list()
  
  for(x in rownames(GFC)){
    
    temp <- as.data.frame(t(GFC[rownames(GFC) %in% c(x),]))
    colnames(temp) <- "clust"
    tmp <- data.frame(chr = as.character(rep(x,ncol(GFC))) , start = rep(0,ncol(GFC)), end = rep(1,ncol(GFC)) , value1 = temp$clust)
    bed4[[x]] <- tmp
    
  }
  
  bed4 <- list.rbind(bed4)
  rownames(bed4) <- 1:nrow(bed4)
  colnames(bed4) <- c("chr", "start", "end", "value1")
  
  
  
  tf_col <- colorRamp2(breaks = c(0,1), colors = c("#afaba5","#e89702"))
  epi_col <- colorRamp2(breaks = c(0,1), colors = c("#afaba5","#237704"))
  f = colorRamp2(breaks = range, colors = c("#053061", "white", "#67001F"))
  
  f2 = factor(unique(cluster_information$try))
  
  GFC_list <- list()
  
  for(x in unique(cluster_information$try)){
    
    tmp <- cluster_information[cluster_information$try == x,"Gene"]
    #tmp <- tmp$Gene
    GFC <- as.data.frame(GFC_all_genes[rownames(GFC_all_genes) %in% tmp,1:(ncol(GFC_all_genes)-1)])
    GFC <- apply(GFC,2,mean)
    GFC_list[[x]] <- GFC
  }
  
  GFC <- as.data.frame(list.rbind(GFC_list))
  
  if(hm_vec == TRUE){
    
    GFC <- GFC[heatmap_vec]
    
  }else{} 
  
  df <- data.frame(x=f2, y=number_TF, z = number_Epi)
  
  df <- df[sort(levels(f2)),]
  
  GFC <- as.data.frame(t(GFC))
  
  GFC <- GFC[levels(df$x)]
  
  colnames(GFC) <- df$y.names
  
  GFC <- as.data.frame(t(GFC))
  
  df$g <- GFC
  
  df$names <- paste0(df$x," ","[", df$y.total,"]")
  
  col <- colorRamp2(breaks = c(min(df$y.total), max(df$y.total)), colors = c("white", "black"))
  
  df$color <- col(df$y.total)
  
  library(circlize)
  circos.clear()
  circos.par(gap.degree = 0,start.degree = 90, clock.wise = TRUE, unit.circle.segments =seg, track.margin = c(0, 0))
  circos.initialize(factors = df$y.names, xlim = c(1,length(df$x)))
  circos.trackPlotRegion(factors=df$y.names, ylim = c(0, 1), bg.col = df$color,bg.border = "grey16", track.height = 0.2)
  
  
  g <- circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), mean(ylim) , sector.name, cex=0.65 , facing = "clockwise", niceFacing = TRUE, col = "#7c0734")
    
  }, bg.border = NA )
  
  circos.genomicTrackPlotRegion(bed4,bg.col= "white", stack = T,panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = f(value[[1]]), 
                       border = "lightgrey", lwd = 0.1, posTransform = posTransform.default, ...)
  }, bg.border = "grey16", track.height = 0.05)
  
  circos.genomicTrackPlotRegion(bed2,bg.col= "white", stack =T,panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = epi_col(value[[1]]), 
                       border = epi_col(value[[1]]), lwd = 0.1, posTransform = posTransform.default, ...)
  }, bg.border = "grey16", track.height = 0.1)
  
  
  
  circos.genomicTrackPlotRegion(bed3,bg.col= "white", stack = T,panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = tf_col(value[[1]]), 
                       border = tf_col(value[[1]]), lwd = 0.1, posTransform = posTransform.default, ...)
  }, bg.border = "grey16", track.height = 0.1)
  
  if(y_comp == T){
    
    GO_enrich <- list()
    for(x in unique(as.character(y_compare$Cluster))){
      
      tmp <- as.vector(y_compare[y_compare$Cluster==x,"Description"])
      immune <- as.numeric(length(tmp[tmp %in% GO_immune_response]))
      metabolic <- as.numeric(length(tmp[tmp %in% GO_metabolic_process]))
      signalling <- as.numeric(length(tmp[tmp %in% GO_signalling]))
      tmp2 <- as.vector(y_compare[y_compare$Cluster==x,"GeneRatio"])
      tmp2 <- as.character(tmp2[1])
      size <- as.numeric(unlist(strsplit(tmp2, split = "/", fixed=T))[2])
      df3 <- data.frame(GO_immune_response = as.numeric(immune),GO_metabolic_process =as.numeric(metabolic),
                        GO_signalling = as.numeric(signalling))
      GO_enrich[[x]] <- df3
      
    }
    enriched_df <- list.rbind(GO_enrich)
    
    f2 = unique(cluster_information$try)
    
    df2 <- data.frame(x=as.character(f2))
    df2$x <- as.character(df2$x)
    
    rest <- as.character(df2[!df2$x %in% rownames(enriched_df),"x"])
    
    matrix <- as.data.frame(matrix(rep(0,(length(rest)*3)),ncol=3))
    rownames(matrix) <- rest
    colnames(matrix) <- colnames(enriched_df)
    
    enriched_df <- rbind(as.data.frame(enriched_df), as.data.frame(matrix))
    
    enriched_df <- as.data.frame(t(enriched_df))
    
    enriched_df <- enriched_df[as.character(df$x)]
    
    enriched_df <- as.data.frame(t(enriched_df))
    
    rownames(enriched_df) <- df$y.names
    
    test_df <- data.frame(factors = rep(rownames(enriched_df),3), x= rep(c(10,20,30),each=nrow(enriched_df)), y=as.vector(unlist(enriched_df)), v= rep(c(1,2,3), each=nrow(enriched_df)))
    
    test_df$z <- apply(test_df,1, function(x){
      
      if((as.numeric(x[3])==as.numeric(0))==T){
        
        "white"
        
      }else{
        
        if(x[4]==1){"#80ad03"}else{if(x[4]==2){"#91012c"}else{"#004089"}}
        
      }
      
    })
    
    test_df$z <- as.character(test_df$z)
    
    test_df$names <- rep(df$names,3)
    
    circos.trackPlotRegion(track.index=5,factors=df$names, ylim = c(0, max(test_df$y)), bg.col = "white",bg.border = "grey16", track.height = 0.1)
    
    circos.trackLines(test_df$names, test_df$x[order(test_df$x)], test_df$y[order(test_df$x)], col = test_df$z, lwd=3, type="h")
  }else{
    
    "GO enrichment was left out" 
    
  }
  
}



plot_genes_on_networks<-function(igraph_object=igraph_object,geneset_name = geneset_name, geneset = geneset, print_to_pdf=FALSE){
  

  # igraph_list<-list()
  # igraph_list[[1]]<-igraph_object
  # layout_df<-BiocGenerics::do.call(c("layout_with_lgl"), igraph_list)
  # 
  node_color <- cluster_information
  node_color$marked_color <- ""
  node_color[node_color$Gene %in% geneset,]$marked_color <- "black"
  node_color[!node_color$Gene %in% geneset,]$marked_color <- rgb(0,0,0,.1)

  if(print_to_pdf==TRUE){
    pdf(paste0("Genes_on_networks_", geneset_name,".pdf"), width = 8 , height = 8)
    
    plot(igraph_object,layout=layout,
         vertex.label=NA,
         edge.color="lightgrey",
         vertex.shape="fcircle",
         vertex.size=3,
         vertex.color=node_color$marked_color,
         vertex.frame.color=node_color$vertex_color,
         vertex.frame.width=0.001)
    
    dev.off()
  }else{
    plot(igraph_object,layout=layout,
         vertex.label=NA,
         edge.color="lightgrey",
         vertex.shape="fcircle",
         vertex.size=3,
         vertex.color=node_color$marked_color,
         vertex.frame.color=node_color$vertex_color,
         vertex.frame.width=0.001)
  }
}














#clusterprofiler
clusterprofiler_GO<-function(cluster_to_check=c("all"),group=c("merged")){
  summary<-allargs()
  cluster_available<-unique(cluster_information$try)
  df<-data.frame()
  no<-0
  biggest_cluster <- cluster_information[grepl(paste0("^" , names(which.max(table(cluster_information$try))) , "$") , cluster_information$try) , ]
  df <- data.frame(longest_cluster = as.character(cluster_information[cluster_information$try == unique(biggest_cluster$try) , 2]))
  
  if(cluster_to_check==c("all")){
    for(Thomas in cluster_available) {
      no <- 1 + no
      x <- as.character(Thomas)
      new.col <- as.character(cluster_information[cluster_information$try == x , 5])
      df[ , no] <- c(new.col , rep(NA , nrow(df)-length(new.col)))
      colnames(df)[no] <- Thomas
    }
  }else{
    for(Thomas in cluster_to_check){
      no <- 1 + no
      x <- as.character(Thomas)
      new.col <- as.character(cluster_information[cluster_information$try == x , 5])
      df[ , no] <- c(new.col , rep(NA , nrow(df)-length(new.col)))
      colnames(df)[no] <- Thomas
    }
  }
  
  ###cluster Profiler
  normdata<- as.data.frame(t(Dataset_1))
  normdata$merged<-info_Dataset[,group]
  universe <- rownames(Dataset_1)
  
  cluster_genes_original<-data.frame(Gene=df)
  returnerlist<-list()
  
  if(!(base::exists(c("cutoff_wd")))){
    cutoff_wd<-file.path(originalwd,chosen_cutoff)
  }
  
  x <- cutoff_wd
  dir.create(x)
  setwd(x)
  plotPath = file.path(x, "clusterProfiler");
  dir.create(file.path(x, "clusterProfiler"), showWarnings = FALSE)
  
  if(organism==c("Mouse")){
    orga_type= c("mouse")
  }else{
    orga_type=c("human")
  }
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("biomaRt")
  #library(biomaRt)
  
  #human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  if(orga_type == "mouse"){
    universe_mouse_human = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = universe, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    universe_mouse_human <- universe_mouse_human[,2]
    universe_Entrez_mouse = clusterProfiler:: bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    universe_Entrez_mouse = unlist(universe_Entrez_mouse[2],use.names = FALSE)
    universe_Entrez_mouse_human = clusterProfiler::bitr(universe_mouse_human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez_mouse_human = unlist(universe_Entrez_mouse_human[2],use.names = FALSE)
  }else{
    universe_Entrez = clusterProfiler:: bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    universe_Entrez = unlist(universe_Entrez[2],use.names = FALSE)
  }
  c1_hallmark_genes <- gmtfile
  list_of_entrez <- list()
  i <- 1
  p <- list()
  #######
  for(id_it in 1:ncol(cluster_genes_original)){
    
    
    
    list_of_genes <- list(cluster_genes_original[,id_it])
    
    #if mouse, mouse symbols are translated to human, because some gene sets works only with human
    if(orga_type == "mouse"){
      universe_orignal <- universe
      genes_cluster = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unlist(list_of_genes), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      cluster_genes <- genes_cluster[,2]
      list_of_genes_mouse_human <- list(cluster_genes)
      entrez_de = clusterProfiler :: bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
      module_entrez_mouse <- unlist(entrez_de[2],use.names = FALSE)
      entrez_de = clusterProfiler::bitr(unlist(list_of_genes_mouse_human), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez_mouse_human <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez_mouse
      
    }else{
      entrez_de = clusterProfiler::bitr(unlist(list_of_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      module_entrez <- unlist(entrez_de[2],use.names = FALSE)
      
      list_of_entrez[[colnames(cluster_genes_original)[id_it]]] <- module_entrez
    }
    
    
    if(orga_type == "mouse"){
      enrichGO <- clusterProfiler::enrichGO(gene = module_entrez_mouse,
                                            universe = universe_Entrez_mouse,
                                            OrgDb = org.Mm.eg.db,
                                            #keytype = 'ENTREZID',
                                            ont = "BP",
                                            pAdjustMethod = "none",
                                            pvalueCutoff  = 0.01,
                                            qvalueCutoff  = 1,
                                            readable      = T)
      
    }else{
      enrichGO <- clusterProfiler::enrichGO(gene = module_entrez,
                                            universe = universe_Entrez,
                                            OrgDb = org.Hs.eg.db,
                                            #keytype = 'ENTREZID',
                                            ont = "BP",
                                            pAdjustMethod = "none",
                                            pvalueCutoff  = 0.05,
                                            qvalueCutoff  = 1,
                                            readable      = T)
      
    }
    
    
  
  
      
      p[[i]] <- barplot(enrichGO, font.size = 24, title = "GO enrichment",showCategory=6)
      
    
      i <- i+1
    
    
  }
  
  do.call(grid.arrange,p)
  do.call("grid.arrange", c(p, ncol=1))
  
}



