
# library(qgraph)
# library(stringr)
# library(magrittr)
# library(dbscan)
# library(STRINGdb)
# library(scales)


# function to extract genes classiefied as hubs
extract_hub_genes <- function(component, top_percentage_for_hubs){

  # extract the degree of each gene ( = vertex) in the component. "degrees" is a data frame with two columns: "degree" and "gene"
  degrees <- data.frame(degree = igraph::degree(component, v = V(component), mode=c("all"), normalized = FALSE)) %>%
    tibble::rownames_to_column("gene")


  # determine the degree quantile that marks the minimum number of degrees in order for a node considered to be a hub based on
  # the parameter "top_percentage_for_hubs" set in the "create_I_Gin" function:
  quant <- quantile(degrees[, "degree"], probs = c(1 - top_percentage_for_hubs), na.rm = TRUE)


  # extract hub genes
  hub_genes <- dplyr::filter(degrees, degree >= quant)[, "gene"]


  return(hub_genes)
}







# function to extract all edges incident to hub genes
extract_hub_edges <- function(edges, hub_genes, allowed_edges, allowed_edges_between_hubs){

  #data frame to store edges incident to hubs
  hub_edges <- NULL

  for (gene in hub_genes){
    #temporary data frame to store edges that contain current hub gene. Later added to "hub_genes".
    tmp <- dplyr::filter(edges, V1 == gene | V2 == gene)

    #detect edges that connect two hubs in order to remove them if their number exceeds the set threshold (allowed_edges_between_hubs)
    hub_to_hub <- dplyr::filter(tmp, V1 %in% hub_genes & V2 %in% hub_genes)

    #remove all hub to hub edges from tmp and then add only the number of allowed ones with highest correlation
    tmp <- dplyr::anti_join(tmp, hub_to_hub, by = c("V1", "V2", "rval")) %>%
      rbind(hub_to_hub[1:allowed_edges_between_hubs,]) %>%
      dplyr::arrange(dplyr::desc(rval)) %>%
      na.omit()

    #reduce the hubs degree to the allowed degree ("allowed_edges")
    if (tmp %>% nrow() > allowed_edges){
      tmp <- tmp[1:allowed_edges,]
    }

    #add tmp to final data frame
    hub_edges <- hub_edges %>% rbind(tmp)

  }

  # remove duplicates
  hub_edges <- dplyr::distinct(hub_edges)


  # add both versions of merged gene names for each edge. This is needed later for identification of known interactions using
  # STRINGdb
  hub_edges$merged <- base::paste0(hub_edges$V1, hub_edges$V2)
  hub_edges$merged_rev <- base::paste0(hub_edges$V2, hub_edges$V1)


  # add aesthetics to the data frame that are needed later for plotting
  hub_edges$color <- "grey"
  hub_edges$curved <- "FALSE"


  return(hub_edges)

}







# function to add node labels
label_nodes <- function(nodes, edges, top_percentage_for_hubs, organism, normal_color, label_TFs = F, color_TF = NA){

  #initialise coloumns needed to set respective aesthetics:
  nodes$label = NA
  nodes$color_label = normal_color


  # label transcription factors
  if(label_TFs){

    if (organism %in% c("human", "Human")){

      # filter genes for transcription factors
      tmp <- dplyr::filter(nodes, gene %in% as.character(TF_list$Human))


      # reset the label and the color of all TFs
      nodes <- dplyr::anti_join(nodes, tmp, by = c("gene", "color", "label", "color_label"))
      tmp$label <- tmp$gene
      tmp$color_label <- color_TF
      nodes <- dplyr::bind_rows(nodes, tmp)

    }else if (organism %in% c("mouse", "Mouse")){

      # filter genes for transcription factors
      tmp <- dplyr::filter(nodes, gene %in% as.character(TF_list$Mouse))

      # reset the label and the color of all TFs
      nodes <- dplyr::anti_join(nodes, tmp, by = c("gene", "color", "label", "color_label"))
      tmp$label <- tmp$gene
      tmp$color_label <- color_TF
      nodes <- dplyr::bind_rows(nodes, tmp)
    }

  }


  # label hub genes
  labels <- table(c(as.character(edges$V1), as.character(edges$V2))) %>%
    BiocGenerics::sort(., decreasing = T) %>%
    as.data.frame(.)

  quant <- quantile(labels[,"Freq"],probs = c(1 - top_percentage_for_hubs))
  labels <- dplyr::filter(labels, Freq >= quant)
  nodes <- base::within(nodes, label[gene %in% as.character(labels$Var1)] <- gene[gene %in% as.character(labels$Var1)])

  # return the modified nodes data frame
  return(nodes)
}







# function to set the edge colours according to known and unknown interactions
set_colours <- function(interactions, edges){

  # iterate over the graph's edges
  for (line in 1:nrow(edges)){

    # if interaction is known, colour the edge accordingly on set the curved aesthetic to FALSE
    if((edges[line, "merged"] %in% interactions$merged) | (edges[line, "merged_rev"] %in% interactions$merged)){

      edges[line, "color"] <- c("#c95555")
      edges[line, "curved"] <- FALSE

      # otherwise set colour to be grey
    }else{

      edges[line, "color"] <- "grey"

    }
  }

  #return the modified edge data frame
  return(edges)
}







# function to translate STRING gene IDs to gene names
translate_gene_names <- function(interactions, mapping){

  interactions$from <- mapping$gene[match(interactions$from, mapping$STRING_id)]
  interactions$to <- mapping$gene[match(interactions$to, mapping$STRING_id)]

  #return updated interactions data frame
  return(interactions)
}







#function to run the string db query:
run_string <- function(genes, edges, string_threshold){ #genes needs to be a data frame with genes in column

  # set up organism specific data base
  if(organism %in% c("mouse", "Mouse")){

    string_db <- STRINGdb$new(version = "10" , species = 10090, score_threshold = string_threshold)

  }else if(organism %in% c("human", "Human")){

    string_db <- STRINGdb$new(version = "11" , species = 9606, score_threshold = string_threshold)

  }else{
    # in case the organism is neither human nor mouse, the ID needs to be entered manually. For the ID, please refer
    # to information on STRING.
    species_ID <- base::readline(promt = "please enter the STRINGdb species ID of your organism: ") %>%
      as.integer()

    string_db <- STRINGdb$new(version = "11" , species = species_ID, score_threshold = string_threshold)

  }

  # status update since this step has a slightly higher time consumption
  print("mapping gene IDs to STRING")

  # mapping the input genes to string identifiers
  mapping <- string_db$map(genes %>% data.frame(gene = ., stringsAsFactors = F), c("gene"), removeUnmappedRows = T)

  # status update since this step has a slightly higher time consumption
  print("")
  print(c("retrieving interactions"))

  # retrieve interactions between input genes from data base
  interactions <- string_db$get_interactions(mapping$STRING_id) %>%
    translate_gene_names(. , mapping = mapping)

  # merge gene names for easier setting of aesthetics in subsequent steps
  interactions$merged = base::paste0(interactions$from, interactions$to)

  # colour edges that represent found interactions
  edges <- set_colours(interactions = interactions, edges = edges)

  #return modified edges data frame
  return(edges)

}







# function to untangle graphs with high density
reduce_edges <- function(edges, hub_genes, max_degree){

  # retrieve non-hub genes
  non_hubs <- c(as.character(edges$V1), as.character(edges$V2)) %>%
    base::unique() %>%
    data.frame(gene = .) %>%
    dplyr::filter(., !gene %in% hub_genes)

  # filter all edges the current non-hub gene is part of, spare those that represent known interactions,
  # sort them by decreasing correlation value
  for (gene in non_hubs$gene){
    tmp <- dplyr::filter(edges, V1 == gene | V2 == gene) %>%
      dplyr::filter(., color == "grey") %>%
      dplyr::arrange(dplyr::desc(rval))

    # only keep the maximal allowed number of edges, chose those with strongest correlation
    edges <- anti_join(edges, tmp, by = c("V1", "V2", "rval", "merged", "merged_rev", "color", "curved")) %>%
      rbind(., tmp[1:max_degree,])

  }

  # return updated edge data frame
  return(edges)
}







# function to set the node size aesthetic:
set_size <- function(node_info, allowed_edges){

  # initialize aesthetic
  node_info$size <- 1

  # scale size of hub nodes with regard to their degree
  node_info <- within(node_info, size[is.na(label) == F] <- as.integer(scales::rescale(degree[is.na(label) == F], to = c(5,20))))

  # return updated data frame
  return(node_info)
}


node_colour <- function(bayes_edges, i_gin_nodes, hub_genes){
  i_gin_nodes$parents <- NA
  i_gin_nodes$children <- NA
  i_gin_nodes$ratio <- NA
  rbPal <- colorRampPalette(c('#FFFFFF','#982D80FF'))(50)[as.numeric(cut(0:100,breaks = 50))]

  for (x in hub_genes){
    i_gin_nodes[i_gin_nodes$gene == x, "parents"] <- dplyr::filter(bayes_edges, to == x) %>%
      nrow()
    i_gin_nodes[i_gin_nodes$gene == x, "children"] <- dplyr::filter(bayes_edges, from == x) %>%
      nrow()
    if(i_gin_nodes[i_gin_nodes$gene == x, "children"]+i_gin_nodes[i_gin_nodes$gene == x, "parents"] == 0){
      i_gin_nodes[i_gin_nodes$gene == x, "ratio"] <- NA
      i_gin_nodes[i_gin_nodes$gene == x, "color"] <- "grey"
    }else{
      i_gin_nodes[i_gin_nodes$gene == x, "ratio"] <- i_gin_nodes[i_gin_nodes$gene == x, "children"] /
        (i_gin_nodes[i_gin_nodes$gene == x, "children"]+i_gin_nodes[i_gin_nodes$gene == x, "parents"])
      #print(rbPal[as.integer(i_gin_nodes[i_gin_nodes$gene == x, "ratio"] *100+1)])
      #print(as.integer(i_gin_nodes[i_gin_nodes$gene == x, "ratio"] *100+1))
      i_gin_nodes[i_gin_nodes$gene == x, "color"] <- rbPal[as.integer(i_gin_nodes[i_gin_nodes$gene == x, "ratio"] *100+1)]
    }

    #print(i_gin_nodes[i_gin_nodes$gene == x, ])
  }

  return(i_gin_nodes)

}




# function to generate and plot the I-Gin for each cluster
plot_components <- function(comps,
                            edges,
                            allowed_edges,
                            allowed_edges_between_hubs,
                            top_percentage_for_hubs,
                            string_threshold,
                            string,
                            max_degree_of_non_hubs,
                            cluster,
                            resampling){


  # collection of igraph objects for all chosen clusters. This will be returned by this function.
  i_gin <- list()


  # iterate over the graph components (corresponding to clusters) in order to create seperate plots:
  for (x in 1:length(comps)){


    #filter for hub genes. "hub_genes" is a vector:
    hub_genes <- extract_hub_genes(component = comps[[x]], top_percentage_for_hubs = top_percentage_for_hubs)


    # extract the edges that are incident to the hub genes
    hub_edges <- extract_hub_edges(edges = edges[[x]], hub_genes = hub_genes,
                                   allowed_edges = allowed_edges, allowed_edges_between_hubs = allowed_edges_between_hubs)


    # extract the gene names of all genes included in the plot (hubs + those connected to them)
    i_gin_genes <- c(as.character(hub_edges$V1), as.character(hub_edges$V2)) %>% base::unique()%>%
      data.frame(stringsAsFactors = F) %>% magrittr::set_colnames("gene")

    # set default colour aesthetic. The nodes will later be re-coloured based on predicted regulatory strenght
    i_gin_genes$color <- "white"


    # add labels and their colours to the data frame, which are needed for the plots aesthetics
    i_gin_genes <- label_nodes(nodes = i_gin_genes, edges = hub_edges, top_percentage_for_hubs = top_percentage_for_hubs,
                               organism = organism, normal_color = c("#231151FF"), label_TFs = T, color_TF = c("#D3436EFF"))


    # if the string parameter was set to TRUE, a STRINGdb query will be conducted:
    if(string){
      hub_edges <- run_string(genes = i_gin_genes$gene, edges = hub_edges, string_threshold = string_threshold)
    }


    # to untangle the graphs structure, non-hub genes can only be addressed by a set maximum of edges from hub genes. Those
    # edges with the highest correlation. The default is 3.
    hub_edges <- reduce_edges(edges = hub_edges, hub_genes = hub_genes, max_degree = max_degree_of_non_hubs) %>% na.omit()

    output <- list()
    bayes_out <- Bayes(hub_genes, hub_edges, resampling = resampling)
    output[["bayes_out"]] <- bayes_out

    i_gin_genes <- node_colour(bayes_edges = bayes_out, i_gin_nodes = i_gin_genes, hub_genes = hub_genes)
    # remove the removed edges from the graph as well
    i_gin[[x]] <- igraph::graph_from_data_frame(hub_edges, vertices = NULL, directed = FALSE )


    # update degrees after edge reduction
    i_gin_genes$degree <- igraph::degree(i_gin[[x]], v = V(i_gin[[x]]), mode=c("all"), normalized = FALSE)


    # set node sizes based on their degrees
    i_gin_genes <- set_size(node_info = i_gin_genes, allowed_edges = allowed_edges)

    windowsFonts("Helvetica" = windowsFont("Helvetica"))
    # set the aesthetics for the igraph object
    igraph_options(vertex.color = i_gin_genes$color,
                   vertex.label = i_gin_genes$label,
                   vertex.label.color = i_gin_genes$color_label,
                   vertex.size = i_gin_genes$size,
                   vertex.label.family="Helvetica",
                   vertex.label.font=2,
                   vertex.label.cex=0.5,
                   vertex.frame.color = i_gin_genes$color_label,
                   edge.color = hub_edges$color,
                   edge.curved = hub_edges$curved)


    # create final graph
    i_gin[[x]] <- igraph::graph_from_data_frame(hub_edges, vertices = NULL, directed = FALSE )


    # calculate the graph layout
    l <- qgraph::qgraph.layout.fruchtermanreingold(igraph::get.edgelist(i_gin[[x]],names = FALSE),vcount=igraph::vcount(i_gin[[x]]),
                                                   weights = igraph::edge.attributes(i_gin[[x]])$weight,
                                                   area=10*(igraph::vcount(i_gin[[x]])^2*2),#2*2
                                                   repulse.rad=(igraph::vcount(i_gin[[x]])^3.1),    #3.1
                                                   niter = 1000)


    # plot I-Gin
    pdf(file=paste0("I_Gin_", cluster[x], ".pdf"))
    #pdf(file=paste0(working_directory, global_settings$save_folder, "/I_Gin_", cluster[x], ".pdf"))
    plot(i_gin[[x]], layout = l, main = cluster[x])
    dev.off()


  }
  output[["i_gin"]] <- i_gin
  # return collection of graph objects for all selected clusters
  return(output)

}







create_I_Gin <- function(data,                  # filtered correlations (above cut-off)
                         cluster_names,         # vector of cluster names (strings) or "all" for all clusters
                         cluster_info,
                         top_percentage_for_hubs = 0.25,
                         allowed_edges = 25,
                         allowed_edges_between_hubs = 3,
                         max_degree_of_non_hubs = 3,
                         string = T,
                         string_threshold = 200,
                         label_all_TF = TRUE,
                         resampling){


  # in case all clusters are to be plotted, extract cluster names except white ("trash clsuter")
  if(cluster_names[1] == "all"){
    cluster_names = cluster_info$color[cluster_info$color != "white"] %>% base::unique()
  }

  edges <- list()
  comps <- list()

  for (c in cluster_names){

    # extract gene names of selected cluster. "genes" is a character vector.
    genes <- dplyr::filter(cluster_info, color == c)[,"gene_n"] %>%
      dplyr::pull() %>%
      stringr::str_split(pattern = ",") %>%
      purrr::flatten_chr()


    # extract all gene pairs where both are in the cluster to be plotted.
    # "edges" is a data frame with the same structure as "data", sorted by r-value in descending manner.
    edges[[c]] <- dplyr::filter(data, V1 %in% genes & V2 %in% genes) %>%
      dplyr::arrange(dplyr::desc(rval))


    #create an i grpah object from edges
    comps[[c]] <- igraph::graph_from_data_frame(edges[[c]], vertices = NULL, directed = FALSE)


  }

  # create an plot an igraph object for each cluster and return the collection of igraph objects
  i_gin <- plot_components(comps = comps, edges = edges, allowed_edges = allowed_edges,
                           allowed_edges_between_hubs = allowed_edges_between_hubs,
                           top_percentage_for_hubs = top_percentage_for_hubs,
                           string_threshold = string_threshold, string = string,
                           max_degree_of_non_hubs = max_degree_of_non_hubs,
                           cluster = cluster_names,
                           resampling = resampling)

  return(i_gin)
}

















