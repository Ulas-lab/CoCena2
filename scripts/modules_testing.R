library("reshape2")

make_rownames_unique <- function(counts_new){

  counts_new <- counts_new[!duplicated(counts_new$SYMBOL),] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., "SYMBOL")
  # remove all non-numeric columns (description, gene id, etc.):
  for (x in colnames(counts_new)){
    if(!is.numeric(counts_new[[x]])){
      counts_new[[x]] <- NULL
    }
  }
  if(ncol(counts_new) == 0){
    print("All columns were deleted when removing non-numeric columns.
          Please check the data type of your expression values.")
  }

  return(counts_new)
  }

# Gene names will be used as rownames, please provide a column names "SYMBOL" !
read_expression_data <- function(file, rown = T, sep = "\t"){
  if(rown){
    expression_data <- read.table(file = file, row.names = 1,
                                  stringsAsFactors = F, sep = sep, check.names = F)
  }else{
    expression_data <- read.table(file = file,
                                  stringsAsFactors = F, sep = sep, check.names = F)
  }

  expression_data <- make_rownames_unique(expression_data)
  return(expression_data)
}


# Sample IDs will be used as rownames, please provide a column names "SampleID" !
read_anno <- function(file, rown = T, sep = "\t"){
  if(rown){
    anno <- read.table(file = file, row.names = 1,
                       stringsAsFactors = F, sep = sep, check.names = T)

  }else{
    anno <- read.table(file = file,
                       stringsAsFactors = F, sep = sep, check.names = T)
  }

  rownames(anno) <- anno$SampleID
  return(anno)
}


calculate_GFCs <- function(expressions, anno, genes, grpvar, range_GFC){ # gene of which GFC shall be calculated

  # filter expression data for genes of which the GFC shall be calculated
  expressions <- tibble::rownames_to_column(expressions, var = "SYMBOL")%>%
    dplyr::filter(., SYMBOL %in% genes)%>%
    tibble::column_to_rownames(., "SYMBOL")

  # calculate GFCs for groups
  if(!grpvar %in% colnames(anno)){
    print("The column you specified to contain the grouping variables does not exists in the annotation table.")
    return()
  }

  #g <- "SMAD5"
  GFCs <- NULL
  for (g in genes){

    if(g %in% rownames(expressions)){
      g_df <- expressions[rownames(expressions) == g,]
    colnames(g_df) <- colnames(expressions)

    # vector of group means:
    mean_vec <- NULL
    for (var in unique(anno[[grpvar]])){
      var_anno <- anno[anno[[grpvar]] == var,]
      var_exp <- g_df[, colnames(g_df) %in% var_anno$SampleID]
      if(is.data.frame(var_exp)){
        mean_vec <- c(mean_vec,mean(var_exp[1,] %>% as.numeric()))
      }else{
        mean_vec <- c(mean_vec,var_exp)

      }

    }
    # overall expression mean of current gene:
    g_mean <- mean(mean_vec)
    GFC_vec <- NULL
    for (x in 1:length(mean_vec)){
      GFC_vec <- c(GFC_vec, gtools::foldchange(mean_vec[x], g_mean) %>%
                     ifelse(.>range_GFC, range_GFC,.) %>%
                     ifelse(.< (-range_GFC), -range_GFC,.))

    }
    # add this genes GFC vector to GFCs
    GFCs <- rbind( GFCs, GFC_vec)
    }
  }

  GFCs <- as.data.frame(GFCs) %>%
    round(., digits = 3)
  # add column with gene names
  GFCs$Gene <- genes
  colnames(GFCs) <- c(unique(anno[[grpvar]]), "Gene")
  rownames(GFCs) <- GFCs$"Gene"
  return(GFCs)

}


cluster_partition_and_plot <- function(cluster_info, GFCs,clust_r,clust_c, order = order){

  # filter for included clusters (non-white)
  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  mat_heatmap <- NULL

  for (c in unique(c_df$color)){
    #get genes from the original cluster
    genes <- c_df[c_df$color == c, ] %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(., split = ",") %>%
      unlist(.)

    genes <- intersect(genes,rownames(sample_file))



    # GFCs of new data set, where genes are found in original cluster
    c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
    c_GFC_means <- apply(c_GFCs[, c(1:(ncol(c_GFCs)-1))], 2, mean)

    mat_heatmap <- rbind(mat_heatmap, c_GFC_means)

  }
  rownames(mat_heatmap) <- c_df$color

  cluster_colors <- factor(c_df$color)
  names(cluster_colors) <- c_df$color





  ha <- ComplexHeatmap::HeatmapAnnotation(clusters = anno_simple(c_df$color, col = cluster_colors,
                                                                 simple_anno_size = unit(0.5, "cm")),
                                          which = "row",
                                          width = unit(0.5, "cm"),
                                          annotation_name_side = "top",
                                          gap = unit(2, "mm"),
                                          annotation_name_rot = 0,
                                          annotation_name_gp = gpar(fontsize = 8))

  hm <- ComplexHeatmap::Heatmap(mat_heatmap,
                                right_annotation = ha,
                                col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1))),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = clust_c,
                                cluster_rows = clust_r,
                                column_names_rot = 90,
                                column_names_centered = F,
                                row_names_gp = gpar(fontsize = 10),
                                rect_gp = grid::gpar(col = "grey"),
                                column_order = order,
                                heatmap_legend_param = list(title = "", legend_height = unit(3, "cm")))

  print(hm)
  return(hm)

}


compare_external_signature <- function(sample_file,
                                       anno_file,
                                       grpvar,
                                       cluster_info,
                                       range_GFC,
                                       clust_r = T,
                                       clust_c = T,
                                       order=order_vec){


  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  genes <- NULL
  for (c in unique(c_df$color)){
    #get genes from the original cluster
    genes <- c(genes, c_df[c_df$color == c, ] %>%
                 dplyr::pull(., "gene_n") %>%
                 base::strsplit(., split = ",") %>%
                 unlist(.))
  }

  genes <- intersect(genes,rownames(sample_file))

  cluster_info_long_other <- as.data.frame(cluster_info_long_format(cluster_information))
  cluster_info_long_other <- cluster_info_long_other[cluster_info_long_other$gene_n %in% genes,]
  cluster_info_long_other <- cluster_info_long_other[cluster_info_long_other$color != "white",1]
  cluster_info_long_other <- as.data.frame(table(cluster_info_long_other))


  cluster_info_long_our <- as.data.frame(cluster_info_long_format(cluster_information))
  cluster_info_long_our <- cluster_info_long_our[cluster_info_long_our$color != "white",1]
  cluster_info_long_our <- as.data.frame(table(cluster_info_long_our))

  cluster_info_long_combined <- merge(cluster_info_long_our,cluster_info_long_other, by.x = "cluster_info_long_our",by.y = "cluster_info_long_other")
  colnames(cluster_info_long_combined) <- c("module","our","other")

  cluster_info_long_combined_m <- melt(cluster_info_long_combined)

  barp <- ggplot( cluster_info_long_combined_m, aes(x=module, y=value, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    coord_flip() +
    xlab("") +
    theme_bw()


  tmp_GFCS <- calculate_GFCs(expressions=sample_file, anno=anno_file, genes =  genes, grpvar = grpvar, range_GFC = range_GFC)
  tmp_mat <- cluster_partition_and_plot(cluster_info, tmp_GFCS,clust_r = clust_r, clust_c = clust_c,order=order)

  return(barp)

}
