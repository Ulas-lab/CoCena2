# data$lnc_counts <- data$lnc_counts[rownames(data$lnc_counts) %in% data$gencode_lnc$gene_name,]
# #install.packages("bnstruct")
# library(bnstruct)
#
# library(readr)


Signordb_human <- read_delim(paste0(working_directory,"scripts/", "Signordb_human.tsv"),
                              "\t", escape_double = FALSE, trim_ws = TRUE)

get_node_expression_vals <- function(edgelist){
  # genes <- dplyr::filter(cluster_calc$cluster_information, color == cluster)%>%
  #   dplyr::pull(., gene_n) %>%
  #   base::strsplit(., split = ",") %>%
  #   BiocGenerics::unlist(.)
  genes <- c(edgelist$V1 %>% as.character(), edgelist$V2 %>% as.character())

  cluster_exp <- subset(count_file_name, rownames(count_file_name) %in% genes) #hier einfach deine gene expression table auf die Zeilen
	#filtern, die den Genen in "genes" entsprechen


  return(cluster_exp)
}

# cluster_exp <- get_cluster_expression_vals("darkgreen")
# bootstrap_input <- list(ahrr = cluster_exp[,1:3], wt = cluster_exp[, 4:6])
# bootstrap_output <- NULL
# for (x in bootstrap_input){
#   tmp <- apply(x,1,sample,
#                size = 8000,
#                replace = T)
#   bootstrap_output <- rbind(bootstrap_output, tmp)
# }
# boot_t <- as.data.frame(bootstrap_output)
# test <- boot_t
#
# d_from_d <- BNDataset(data = test,
#                       discreteness = rep(F, ncol(test)),
#                       variables = colnames(test),
#                       node.sizes = rep(6, ncol(test)))
#
# #
# # net1 <- learn.network(d_from_d,
# #                       algo = "mmhc",
# #                       scoring.func = "BIC",
# #                       bootstrap = F,
# #                       max.fanin = 1,
# #                       use.cpc = T,
# #                       use.imputed.data = F)
# #
# # plot(net1)
#
# known_causalities <- dplyr::filter(Signordb_human, ENTITYA %in% base::toupper(rownames(cluster_exp)) &
#                                      ENTITYB %in% base::toupper(rownames(cluster_exp)) & !ENTITYA == ENTITYB)
#
# known_g <- igraph::graph_from_edgelist(cbind(known_causalities$ENTITYA, known_causalities$ENTITYB), directed = T)
# known_g <- as_adjacency_matrix(known_g, names = TRUE, sparse = FALSE, attr = NULL)
# known_g[!known_g == 0] <- 1
# cluster_known <- t(test) %>% as.data.frame()
# cluster_known <- cluster_known[base::toupper(rownames(cluster_known)) %in% rownames(known_g),]
# cluster_known <- t(cluster_known) %>% as.data.frame()
# d_from_known <- BNDataset(data = cluster_known,
#                       discreteness = rep(F, ncol(cluster_known)),
#                       variables = colnames(cluster_known),
#                       node.sizes = rep(6, ncol(cluster_known)))
# testbn <- bnstruct::BN(data = d_from_known)
# testbn@dag <- known_g
#
# mandatory <- matrix(rep(0, ncol(test)*ncol(test)), ncol = ncol(test))
# rownames(mandatory) <- base::toupper(colnames(test))
# colnames(mandatory) <- base::toupper(colnames(test))
#
# indices<- data.frame(row = row(known_g)[which(!known_g == 0)], col = col(known_g)[which(!known_g == 0)])
# indices$rowname <- NA
# indices$colname <- NA
#
# for (x in 1:nrow(indices)){
#   indices[x,"rowname"] <- rownames(known_g)[indices[x,1]]
#   indices[x,"colname"] <- rownames(known_g)[indices[x,2]]
#   mandatory[rownames(known_g)[indices[x,1]], rownames(known_g)[indices[x,2]]] <- 1
# }
#
# rownames(mandatory) <- colnames(test)
# colnames(mandatory) <- colnames(test)
#
# net1 <- learn.network(d_from_d,
#                       algo = "mmhc",
#                       scoring.func = "BIC",
#                       bootstrap = F,
#                       max.fanin = 1,
#                       initial.network = mandatory,
#                       use.cpc = F,
#                       use.imputed.data = F)





Bayes <- function(hub_genes, edgelist, resampling){
  bayes_out <- NULL
  for( x in hub_genes){
    tmp_edges <- dplyr::filter(edgelist, as.character(V1) == x | as.character(V2) == x)
    cluster_exp <- get_node_expression_vals(tmp_edges)
    #print(cluster_exp)

    if(resampling == T){
      bootstrap_input <- list(ahrr = cluster_exp[,1:3], wt = cluster_exp[, 4:6])
      bootstrap_output <- NULL
      for (y in bootstrap_input){
        #print(y)
        tmp <- apply(y,1,sample,
                     size = 30,
                     replace = T)
        bootstrap_output <- rbind(bootstrap_output, tmp)
        tmp <- NULL
      }
      boot_t <- as.data.frame(bootstrap_output)
      cluster_exp <- boot_t
    }else{
      cluster_exp <- t(cluster_exp) %>% as.data.frame()
    }

    #print(cluster_exp)
    # cluster_exp <- t(cluster_exp) %>%
    #   as.data.frame()
    d_from_d <- BNDataset(data = cluster_exp,
                          discreteness = rep(F, ncol(cluster_exp)),
                          variables = colnames(cluster_exp),
                          node.sizes = rep(6, ncol(cluster_exp))) #6 needs to be reconsidered

    known_causalities <- dplyr::filter(Signordb_human, ENTITYA %in% base::toupper(colnames(cluster_exp)) &
                                         ENTITYB %in% base::toupper(colnames(cluster_exp)) & !ENTITYA == ENTITYB)
    #print(known_causalities)
    if(nrow(known_causalities) > 0){
      known_g <- igraph::graph_from_edgelist(cbind(known_causalities$ENTITYA, known_causalities$ENTITYB), directed = T)
      known_g <- as_adjacency_matrix(known_g, names = TRUE, sparse = FALSE, attr = NULL)
      known_g[!known_g == 0] <- 1
      cluster_known <- t(cluster_exp) %>% as.data.frame()
      cluster_known <- cluster_known[base::toupper(rownames(cluster_known)) %in% rownames(known_g),]
      cluster_known <- t(cluster_known) %>% as.data.frame()
      #print(cluster_known)
      d_from_known <- BNDataset(data = cluster_known,
                                discreteness = rep(F, ncol(cluster_known)),
                                variables = colnames(cluster_known),
                                node.sizes = rep(6, ncol(cluster_known)))
      testbn <- bnstruct::BN(data = d_from_known)
      testbn@dag <- known_g

      mandatory <- matrix(rep(0, ncol(cluster_exp)*ncol(cluster_exp)), ncol = ncol(cluster_exp))
      rownames(mandatory) <- base::toupper(colnames(cluster_exp))
      colnames(mandatory) <- base::toupper(colnames(cluster_exp))

      indices<- data.frame(row = row(known_g)[which(!known_g == 0)], col = col(known_g)[which(!known_g == 0)])
      indices$rowname <- NA
      indices$colname <- NA

      for (x in 1:nrow(indices)){
        indices[x,"rowname"] <- rownames(known_g)[indices[x,1]]
        indices[x,"colname"] <- rownames(known_g)[indices[x,2]]
        mandatory[rownames(known_g)[indices[x,1]], rownames(known_g)[indices[x,2]]] <- 1
      }

      rownames(mandatory) <- colnames(cluster_exp)
      colnames(mandatory) <- colnames(cluster_exp)

      net1 <- learn.network(d_from_d,
                            algo = "mmhc",
                            scoring.func = "BIC",
                            bootstrap = F,
                            #max.fanin = 2,
                            initial.network = mandatory,
                            #use.cpc = F,
                            use.imputed.data = F)
    }else{
      net1 <- learn.network(d_from_d,
                            algo = "mmhc",
                            scoring.func = "BIC",
                            bootstrap = F,
                            #max.fanin = 2,
                            #use.cpc = F,
                            use.imputed.data = F)
    }

    m <- net1@dag
    rownames(m) <- colnames(m) <- net1@variables
    g <- igraph::graph.adjacency(m)


    bayes_out <- rbind(bayes_out, igraph::get.edgelist(g))
  }
  bayes_out <- as.data.frame(bayes_out)
  colnames(bayes_out) <- c("from", "to")
  return(bayes_out)
}
#
# cluster_exp <- get_node_expression_vals(hub_edges[1:5,])
# print(nrow(cluster_exp))
#
# bootstrap_input <- list(ahrr = cluster_exp[,1:3], wt = cluster_exp[, 4:6])
# bootstrap_output <- NULL
# for (x in bootstrap_input){
#   tmp <- apply(x,1,sample,
#                size = 100,
#                replace = T)
#   bootstrap_output <- rbind(bootstrap_output, tmp)
# }
# boot_t <- as.data.frame(bootstrap_output)
# test <- boot_t
#
# d_from_d <- BNDataset(data = test,
#                       discreteness = rep(F, ncol(test)),
#                       variables = colnames(test),
#                       node.sizes = rep(6, ncol(test)))
#
#
#
#
# plot(net1)
