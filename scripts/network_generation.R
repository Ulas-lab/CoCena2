# Filter igraph object and apply layout algorithm  ----------------------------------------------------------

network_layout <- function (igraph.object,
                            use.layout = c("layout_with_fr", "cytoscape"),
                            min.nodes.number,
                            select.subnetwork = NULL){

  # Generate layout based on full network or only specified clusters

  if(!is.null(select.subnetwork)) {

    vertices_subnetwork <- vertex_meta[vertex_meta$cluster == select.subnetwork, "gene_n"]

    ig <- igraph::induced_subgraph(graph = igraph.object, vids = vertices_subnetwork,
                                   impl = "copy_and_delete")


  } else {

    ig <- igraph.object

    }


  # Generate layout using layout_with_fr or cytoscape (prefuse force-directed)

  num_components = igraph::count_components(ig)

  if(use.layout == "layout_with_fr") {

    if(num_components==1){

      layout_fr <- igraph::layout_with_fr(graph = ig, niter = 500, grid = "grid")
      network_graph <- ig


    } else {
      components_ig <- 1:num_components
      components_numbers <- components_ig[igraph::components(ig)[[2]] >= min.nodes.number] #>
      V(ig)$comp <- igraph::components(ig)[[1]]

      list_of_components_graph <- list()
      list_of_components_graph <- lapply(components_numbers,
                                         function(x) induced_subgraph(ig , V(ig)$comp == x))

      layouts_on_list <- lapply(list_of_components_graph,
                                function(x) layout_with_fr(graph=x,niter=500, grid="grid" )) #LISA: niter/grid

      layout_fr <- merge_coords(graphs = list_of_components_graph , layouts = layouts_on_list)

      network_graph <- disjoint_union(list_of_components_graph)


      #Creates separate graph/layout for components with equal or more than min.nodes.number
      #Merges the graph layouts to get one layout for all components
      #disjoint_union merges the graph objects to obtain one graph for all components


    }

    return_network <- list()
    return_network[[c("graph_object")]] <- network_graph
    return_network[[c("layout")]] <- layout_fr

    return(return_network)

  }


  if(use.layout == "cytoscape") {

    if(num_components == 1) {

      network_graph <- ig


    } else {

      components_ig <- 1:num_components
      components_numbers <- components_ig[igraph::components(ig)[[2]] >= min.nodes.number]
      V(ig)$comp <- igraph::components(ig)[[1]]

      list_of_components_graph <- list()
      list_of_components_graph <- lapply(components_numbers,
                                          function(x)induced_subgraph(ig , V(ig)$comp == x))

      network_graph <- disjoint_union(list_of_components_graph)

    }

    #Modifications to skip time-consuming step (does not work at the moment)
    #source(paste0(working_directory,"scripts/", "RCy3_commands_mod.R" ))

    cytoscapePing()

    #transform igraph to network in cytoscape
    RCy3::createNetworkFromIgraph(network_graph, "MyNetwork")
   # RCy3::layoutNetwork("force-directed")

    #get layout coordinates from cytoscape and sort by igraph vertices
    layout_cytoscape <- RCy3::getNodePosition()
    #
    #  layout_cytoscape$x_location <- as.numeric(levels(layout_cytoscape$x_location)) [layout_cytoscape$x_location]
    #  layout_cytoscape$y_location <- as.numeric(levels(layout_cytoscape$y_location)) [layout_cytoscape$y_location]

    layout_cytoscape$x_location <- as.numeric(layout_cytoscape$x_location)
    layout_cytoscape$y_location <- as.numeric(layout_cytoscape$y_location)
    #

    layout_cytoscape <- layout_cytoscape[names(V(network_graph)),]


    #layout_cytoscape <- as.matrix(layout_cytoscape[names(V(network_graph)),])

    # layout_cytoscape$x <- layout_cytoscape$x_location
    # layout_cytoscape$y <- layout_cytoscape$y_location
#
#     layout_cytoscape$x_location = NULL
#     layout_cytoscape$y_location = NULL
#
     layout_cytoscape <- as.matrix(layout_cytoscape)

    colnames(layout_cytoscape) <- c("x","y")
    #use_layout <- layout_cytoscape

    # ig2 <- createIgraphFromNetwork("MyNetwork")

    #output
    return_network <- list()
    return_network[[c("graph_object")]] <- network_graph
    return_network[[c("layout")]] <- layout_cytoscape

    #igraph_object_pers <- network_graph
    #layout_matrix_pers <- layout_cytoscape

    return(return_network)


  }

}



# Collect node attribute information --------------------------------------------------------------


node_information <- function (igraph.object, data_df, GFC_df, TF_df, hallmark_df, go_df, kegg_df, reactome_df,org) {


  # Preparation of cluster data

  vertex_df = V(igraph.object) %>% names() %>% as.data.frame()
  colnames(vertex_df) ="gene_n"
  assign_size_color_vertex = function(nr, cluster_df){
    dx = cluster_df[nr,]
    df_out=data.frame(gene_n = stri_split_regex(dx$gene_n, pattern = ",")%>% unlist() %>% as.character(),
                      size=dx$vertexsize,
                      color=dx$color,
                      cluster=dx$clusters)
    return(df_out)
  }

  vertex_meta = dplyr::left_join(vertex_df,do.call("rbind",
                                                   lapply(1:nrow(data_df),
                                                          assign_size_color_vertex,
                                                          cluster_df=data_df)))

  colnames(vertex_meta) <- c("Gene", "node_size", "cluster_color", "cluster_name")




  # Preparation of additional info

  TF_df <- TF_df[, c(Hmisc::capitalize(org), "Merged Taxa")]
  colnames(TF_df) <- c("Gene", "TF_info")


  hallmark_collapsed <- hallmark_df %>% group_by_at(vars(gene)) %>%
    summarize_all(paste, collapse=",")
  colnames(hallmark_collapsed) <- c("Gene", "HALLMARK")

  go_collapsed <- go_df %>% group_by_at(vars(gene)) %>%
    summarize_all(paste, collapse=",")
  colnames(go_collapsed) <- c("Gene", "GO")

  kegg_collapsed <- kegg_df %>% group_by_at(vars(gene)) %>%
    summarize_all(paste, collapse=",")
  colnames(kegg_collapsed) <- c("Gene", "KEGG")

  reactome_collapsed <- reactome_df %>% group_by_at(vars(gene)) %>%
    summarize_all(paste, collapse=",")
  colnames(reactome_collapsed) <- c("Gene", "REACTOME")


  GSEA_df <- list(hallmark_collapsed, go_collapsed, kegg_collapsed, reactome_collapsed) %>%
    purrr::reduce(full_join, by = "Gene")



  #Joining
  node_attributes <- list(vertex_meta, GFC_df, TF_df, GSEA_df) %>%
    purrr::reduce(left_join, by = "Gene") %>%
    as.data.frame()


  #Adjustments
  node_attributes$TF_info <- as.character(node_attributes$TF_info)
  node_attributes[is.na(node_attributes$TF_info),"TF_info"] <- "noTF"

  node_attributes$vertex.names <- node_attributes$Gene  #for joining with network object

  return(node_attributes)

}




# Generate Network object for ggplot --------------------------------------------------------------

generate_network_object <- function(graph_object, use_layout) {

  n <- ggnetwork(graph_object, layout = use_layout)
  ggnet_object <- merge(x = n, y = node_attributes, by.x = "name", by.y = "vertex.names") #alternatively: by = c(*insert colname of genes in "x"* = *insert colname of genes in "y"*)

  return(ggnet_object)
}


