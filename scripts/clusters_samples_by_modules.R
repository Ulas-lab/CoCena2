library(ggalluvial)
library(factoextra)

heatmap_clusters_table <- function(data) {

  # Data preparation

   mat_heatmap= data %>%
    dplyr::filter(cluster_included=="yes")%>%
    dplyr::select(grp_means) %>%
    tidyr::separate(col=grp_means, sep=",",
                    convert = T,
                    into = 1:(stri_count_regex(data$grp_means[1], pattern = ",")+1) %>%
                      as.character()) %>% as.matrix()

  rownames(mat_heatmap) =data %>%
    dplyr::filter(cluster_included=="yes") %>%
    dplyr::mutate(labs=color) %>% dplyr::pull(labs)

  colnames(mat_heatmap) =  data %>%
    dplyr::filter(cluster_included=="yes") %>% dplyr::select(conditions) %>% purrr::map(1) %>%
    stri_split_regex(pattern = "#") %>% unlist()


  return(t(mat_heatmap))

}




heatmap_clusters_plot <- function(heatmap_data) {


  plot(hclust(dist(cluster_table)), hang = -1, cex = 0.6)


}

investigate_cluster <- function(cluster_table_var = cluster_table) {

  heatmap_clusters_plot(cluster_table_var)

 # set.seed(1)
  wss <- (nrow(cluster_table_var)-1)*sum(apply(cluster_table_var,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(cluster_table_var,
                                       centers=i)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")

}


heatmap_clusters_set_height <- function(heatmap_data, height) {

  # library(dplyr)
  # library(tidyverse)
  # library(dendextend)
  # library(colormap)
  # library(kableExtra)

  {plot(hclust(dist(cluster_table)), hang = -1, cex = 0.6)
  abline(h=height,col="red")}

  clusters <- as.dendrogram(hclust(dist(cluster_table)))
  new_cluster <- dendextend:::cutree.dendrogram(clusters,h=height)
  k <- as.numeric(sub('.*:', '', summary(new_cluster)[6]))


  fviz_dend(clusters,
            k = k,
            cex = 0.6,                     # Label size
            palette = "jco",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "jco",           # Rectangle color
            labels_track_height = 4      # Augment the room for labels
  )+
    theme(panel.background = element_blank(),
    axis.line.y = element_line(colour = "black", size = 0.5))


}



heatmap_clusters_set_k <- function(heatmap_data = cluster_table, k) {

  clusters <- as.dendrogram(hclust(dist(cluster_table)))

  fviz_dend(clusters,
            k = k,
            cex = 0.6,                     # Label size
            palette = "jco",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "jco",           # Rectangle color
            labels_track_height = 4      # Augment the room for labels
  )+
    theme(panel.background = element_blank(),
          axis.line.y = element_line(colour = "black", size = 0.5))


}

heatmap_insert_new_clustering <- function(heatmap_data, height = NULL, k= NULL, voi_id, compare_column_left,middle,compare_column_right,color) {

  if("new_cluster" %in% colnames(sample_table)){
    sample_table$new_cluster = NULL
    }

  clusters <- as.dendrogram(hclust(dist(cluster_table)))
  new_cluster <- as.data.frame(dendextend:::cutree.dendrogram(clusters,h=height, k=k))
  rownames(new_cluster) <- gsub("GFC_", "", rownames(new_cluster))
  colnames(new_cluster) <- "new_cluster"
  sample_table <- merge(x=sample_table, y=new_cluster, by.x = voi_id , by.y = 0)
  sample_table$new_cluster <- as.factor(sample_table$new_cluster)



  ggplot(as.data.frame(sample_table),
         aes_string(axis1 = compare_column_left, axis2 = middle ,axis3 = compare_column_right)) +
    geom_alluvium(aes_string(fill = color), width = 1/12, color = "black") +
    geom_stratum(width = 1/12, fill = "grey", color = "black") +
    geom_text(stat = "stratum", infer.label = TRUE, reverse = T) +
    scale_x_discrete(limits = c(compare_column_left, middle,compare_column_right), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1")

}
#
# ggplot(as.data.frame(Titanic),
#        aes(y = Freq,
#            axis1 = Survived, axis2 = Sex, axis3 = Class)) +
#   geom_alluvium(aes(fill = Class),
#                 width = 0, knot.pos = 0, reverse = FALSE) +
#   guides(fill = FALSE) +
#   geom_stratum(width = 1/8, reverse = FALSE) +
#   geom_text(stat = "stratum", infer.label = TRUE, reverse = FALSE) +

sample_table_insert_new_clustering <- function(heatmap_data, height, voi_id, table) {

  if("new_cluster" %in% colnames(table)){
    table$new_cluster = NULL
  }

  clusters <- as.dendrogram(hclust(dist(heatmap_data)))
  new_cluster <- as.data.frame(dendextend:::cutree.dendrogram(clusters,h=height))
  rownames(new_cluster) <- gsub("GFC_", "", rownames(new_cluster))
  colnames(new_cluster) <- "new_cluster"
  table <- merge(x=table, y=new_cluster, by.x = voi_id , by.y = 0)
  table$"new_cluster" <- as.factor(table$"new_cluster")
  rownames(table) <- table$ID
  return(table)
}



gene_expression_over_modules <- function(cluster_data=cluster_information,
                                         data_to_test=count_file_name,
                                         genes_of_int=genes_of_interest,
                                         anno_tab = sample_table,
                                         column_of_int = column_of_interest,
                                         order=order_vec,
                                         gene_set_name=NULL,
                                         top_n_modules = NULL){


  cluster_information_l <- cluster_info_long_format(cluster_information)

  if(genes_of_int=="all"){

    count_funct_term <- subset(data_to_test, rownames(data_to_test) %in% cluster_information_l$gene_n)

  }else{
    count_funct_term <- subset(data_to_test, rownames(data_to_test) %in% genes_of_interest)
  }

  count_funct_term_cluster <- merge(count_funct_term,cluster_information_l, by.x = 0, by.y = "gene_n")

  count_funct_term_cluster_m <- melt(count_funct_term_cluster)


  count_funct_term_cluster_m <- merge(count_funct_term_cluster_m,anno_tab[,c("ID",column_of_int)], by.x = "variable", by.y = "ID")


  count_funct_term_cluster_m %>%
    group_by(!!!syms(column_of_int), color,Row.names) %>%
    summarise_at(vars(value), funs(mean(., na.rm=TRUE))) -> count_funct_term_cluster


  if(is.null(top_n_modules)){

    count_funct_term_cluster %>%
      count(color) %>%
      select(color,n) %>%
      group_by(color) %>%
      summarise_at(vars(n), funs(mean(., na.rm=TRUE)))%>%
      mutate(ratio = n/max(n)) -> count_funct_term_cluster_ratio

  }else{

    count_funct_term_cluster %>%
      count(color) %>%
      select(color,n) %>%
      group_by(color) %>%
      summarise_at(vars(n), funs(mean(., na.rm=TRUE)))%>%
      mutate(ratio = n/max(n)) %>%
      top_n(top_n_modules)-> count_funct_term_cluster_ratio

  }

  count_funct_term_cluster <- merge(count_funct_term_cluster, count_funct_term_cluster_ratio, by ="color",all.y = T)

  if(!is.null(order)){
    count_funct_term_cluster[,column_of_int] <- factor(count_funct_term_cluster[,column_of_int], levels = order)
  }

  count_funct_term_cluster <- subset(count_funct_term_cluster,!count_funct_term_cluster$color %in% "white")

  p <- ggplot(count_funct_term_cluster, aes_string(x = column_of_int,  y = "value",color ="color")) +
    geom_line(stat='smooth', aes(fill = color, alpha=ratio, size=ratio,group = color))+
    scale_colour_manual(values = unique(count_funct_term_cluster$color))


  smooth_dat <- setDT(ggplot_build(p)$data[[1]])
  smooth_lab <- smooth_dat[smooth_dat[, .I[x == max(x)], by=group]$V1]

  count_funct_term_cluster$Module <- count_funct_term_cluster$color
  count_funct_term_cluster$Ratio <- count_funct_term_cluster$ratio

  if(is.null(gene_set_name)){
    title_content <- "Gene set expression along modules"
    }else{
      title_content <- paste0(gene_set_name," expression along modules")
      }

  return(count_funct_term_cluster)
  # p <- ggplot(count_funct_term_cluster, aes(x = new_cluster,  y = value,color =Module)) +
  #   geom_line(stat='smooth', aes(fill = Module, size=Ratio,group = Module))+
  #   scale_colour_manual(values = unique(count_funct_term_cluster$Module))+theme_bw()+
  #   geom_text_repel(data = smooth_lab, aes(label = colour, x=x, y=y),inherit.aes = FALSE)+
  #   labs(title= title_content,y="mean expression")
  #
  # return(p)
}



library(fpc)
cstats.table <- function(dist, tree, k) {
  clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                    "wb.ratio","dunn2","avg.silwidth")
  clust.size <- c("cluster.size")
  stats.names <- c()
  row.clust <- c()
  output.stats <- matrix(ncol = k, nrow = length(clust.assess))
  cluster.sizes <- matrix(ncol = k, nrow = k)
  for(i in c(1:k)){
    row.clust[i] <- paste("Cluster-", i, " size")
  }
  for(i in c(2:k)){
    stats.names[i] <- paste("Test", i-1)
    for(j in seq_along(clust.assess)){
      output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
    }
    for(d in 1:k) {
      cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
      dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
      cluster.sizes[d, i]
    }
  }
  output.stats.df <- data.frame(output.stats)
  cluster.sizes <- data.frame(cluster.sizes)
  cluster.sizes[is.na(cluster.sizes)] <- 0
  rows.all <- c(clust.assess, row.clust)
  # rownames(output.stats.df) <- clust.assess
  output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
  colnames(output) <- stats.names[2:k]
  rownames(output) <- rows.all
  is.num <- sapply(output, is.numeric)
  output[is.num] <- lapply(output[is.num], round, 2)
  output
}

