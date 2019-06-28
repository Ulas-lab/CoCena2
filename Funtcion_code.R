#####Set working directory
setwd("E:/Paper/Classifier/2019/")
mainDir <- getwd()
subDir <- c("CoCena")
dir.create(file.path(mainDir , subDir))
setwd(file.path(mainDir , subDir))
originalwd <- getwd()

# load packages

packages_to_load <- c("igraph" , "plotly" ,"ggplot2" ,"bnstruct", "gridExtra" , "Hmisc", "ReactomePA","reactome.db", 
                      "RColorBrewer" , "gtools" , "rlist" , "reshape2" , "gplots" , "moduleColor" , 
                      "NMF" , "rlist" , "ggpubr","RJSONIO","httr","stringr","XML","qgraph","rjson","DOSE","clusterProfiler",
                      "devtools","clues","bnstruct","fields","curl","httpuv","intergraph","bnlearn","viridis","foreign")

lapply(packages_to_load , require , character.only=TRUE)
 #install_github('cytoscape/cytoscape-automation/for-scripters/R/r2cytoscape')
#devtools::install_github("collectivemedia/tictoc")

library(r2cytoscape)
library(igraph)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(gtools)
library(rlist)
library(reshape2)
library(gplots)
library(moduleColor)
library(NMF)
library(rlist)
library(ggpubr)
library(tictoc)
library(grid)
library(gridExtra)
library(viridis)
library(clues)
library(foreign)
library(httpuv)
library(devtools)
library(r2excel)
library(grid)
library(clusterProfiler)
library(ReactomePA)
library(pcaGoPromoter)
library(biomaRt)
library(reshape2)
library(dplyr)

source("E:/ownCloud/backup/R/Cocena_functions.R")

#####################################################################################
#set the organism your data is from
organism=c("Human") #//c("Human")
####################################################################################
library(clusterProfiler)
TF_list <- read.delim(paste0(mainDir,"/TFcat.txt") , header = TRUE , check.names = FALSE)
gmtfile <- clusterProfiler::read.gmt(paste0(mainDir,"/h.all.v6.1.symbols.gmt"))
####################################################################################

tic("Load in data")

#####load in your data

###############################################################################################
#make sure that in your dataset, 
#gene names are as rownames and samples are as column names!
###############################################################################################
Dataset_1 <- read.table("E:/Paper/Classifier/2019/Dataset_1.csv" , header = TRUE , check.names = FALSE, sep = ",", row.names = 1)
#Dataset_1 <- read.table("E:/owncloud/backup/Projects/Viemann_Influenza/kallisto/SVA3/Dataset_1.csv" , header = TRUE , check.names = FALSE, sep = ",", row.names = 1)

###############################################################################################
#make sure that sample names are in rownames
#MAKE sure that there is no " - "in merged names - cytoscape cannot handle these
#Make sure that there is no space in your names- cytoscape cannot handle these
######################################################################################
info_Dataset <- read.table("E:/Paper/Classifier/2019/info_Dataset.csv" , header = TRUE , check.names = FALSE, sep = ",", row.names = 1)
#info_Dataset <- read.table("E:/owncloud/backup/Projects/Viemann_Influenza/kallisto/SVA3/info_Dataset.csv" , header = TRUE , check.names = FALSE, sep = ",", row.names = 1)

info_Dataset <- info_Dataset[,c("Dataset","Condition", "Disease")]
colnames(info_Dataset) <- c("dataset","condition","merged")


####################################################################################



#Save original data to new variable which could be changed during process
#original Dataset can always be found under Dataset_1
original_data <- Dataset_1


#####################################################################################
##filter data if wanted
filter=TRUE
#####################################################################################
#1) on transcription factors
#original_data <- original_data[rownames(original_data) %in% TF_list$Human , ] ## HUman or mouse
#2) any other gene_list

#write.csv(union_DE_genes, "union_DE_genes.csv")
#gene_list<-read.table("E:/ownCloud/backup/Projects/Viemann_Influenza/kallisto/SVA3/union_DE_genes.csv" , header = TRUE , check.names = FALSE, sep = ",", row.names = 1)
#original_data <- original_data[rownames(original_data) %in% gene_list$x , ] ## HUman or mouse
#3) topX
dd1 = original_data[order(apply(original_data, 1, var),decreasing = T),]
original_data <- head(dd1,1000)


toc()

####################################################################################
####################################################################################

tic("Input data heatmap")
collectGarbage()

#####
#first visualisation of your data
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))


pheatmap::pheatmap(original_data, 
                   cluster_row = T,
                   cluster_cols = T,
                   color = col.pal,
                   scale = c("row"),
                   annotation_col = info_Dataset,
                   main = c("Heatmap: Gene Expression Data - whole Dataset "),
                   filename = "Input_data.pdf")

####################################################################################
#There will be a function collecting all set variables values in this list "summary"
###
summary<-list()

toc()



####################################################################################
####################################################################################

tic("Correlation analysis")

###Start CoCena²
correlation_df<-correlation(type_of_correlation="pearson",            #pearson / spearman
                            pVal_treshold=0.05,                       #pValue Treshold for calculated correlation
                            save_raw_data=FALSE)                      #save correlation data in .txt

summary[["correlation"]]<-correlation_df$summary
correlation_df<-correlation_df$raw_data


####deciding on cutoff
library(igraph)
cutoff_options<-cutoff_visualisation(correlation_df=correlation_df,
                                     min_corr=0.5,                    # min correlation 
                                     range_cutoff_length=300,          #No.of cutoffs tested (all higher than min_cor)
                                    min_no_for_cluster=10)            #how many nodes are minimum for a network 

summary[["cutoff_visualisation"]]<-cutoff_options$summary


####################################################################################
##Look at plot
#cutoff_options$Plotly_object
#slide through cutoffs and decide on one
#noch händisch einzutragen
####################################################################################
cutoff_options$Cutoff_df
chosen_cutoff <- optimal_cutoff(cutoff_options)

toc()

####################################################################################

tic("Construct network")
cutoff_wd<-paste0(originalwd,"/",chosen_cutoff)
dir.create(file.path(cutoff_wd))
summary[["chosen_cutoff"]]<-chosen_cutoff
setwd(cutoff_wd)


##plot Network optional change layout
layout_options <- grep("^layout_" , ls("package:igraph") , value = TRUE)[-1]
layout_options
# [1] "layout_as_bipartite"  "layout_as_star"       "layout_as_tree"      
# [4] "layout_components"    "layout_in_circle"     "layout_nicely"       
# [7] "layout_on_grid"       "layout_on_sphere"     "layout_randomly"     
# [10] "layout_with_dh"       "layout_with_drl"      "layout_with_fr"      
# [13] "layout_with_gem"      "layout_with_graphopt" "layout_with_kk"      
# [16] "layout_with_lgl"      "layout_with_mds"      "layout_with_sugiyama"

#####################################################################################
#####test different layouts

layouts_on_list_coord<-test_layout(layouts_to_test=c("layout_with_fr",
                                                     "layout_with_kk",
                                                     "layout_with_lgl"),            #choose from layout_options
                                   min_nodes_number_for_network=10)

summary[["test_layout"]]<-layouts_on_list_coord[[c("summary")]]
layouts_on_list_coord<-layouts_on_list_coord$layouts_on_list_coord
dev.off()


####must run these command to get an igraph object

return_list<-plot_network(data=correlation_df,
                            layout=c("layout_with_kk"),                          #choose anything       if you tested layouts change to layout_on_grid to fasten things up
                            min_nodes_number_for_network=20,
                            show_HC= TRUE)                                       #see HC after cutting whole Dataset

summary[["plot_network"]]<-return_list$summary
igraph_object<-return_list[["graph_object"]]

toc()

######Scale free Topology check
tic("Check scale free topology")

fit_power_law(igraph_object)

####layout can be chosen from
layout<-return_list[["layout"]]

####################################################################################
######if test_layout was run
names(layouts_on_list_coord)# names of your layout options
layout<-layouts_on_list_coord$layout_with_kk
####################################################################################

toc()
#GFC calculation
tic("GFC calculation and visualization")

GFC_all_genes<-GFC_calculation(normdata=original_data,
                               group=c("merged"),       #column name from info_Dataset
                               data_in_log=F,          #is the inloaded data at some point put in log (e.g. DeSeq's rlog)
                               range_GFC=2.0)              #if calculated GFC value is over range_GFC it is 
                                                           #set to range_GFC value due to nice visulaisation
summary[["GFC_calculation"]]<-GFC_all_genes$summary
GFC_all_genes<-GFC_all_genes$GFC_all_genes

plot_GFC_networks(igraph_object=igraph_object,
                  GFC_all_genes = GFC_all_genes,
                  print_to_pdf=T,                      #GFC plots will be plotted in pdf - ATTENTION when further work will be done in Corel Draw not recommended
                  print_edges_png_nodes_pdf=c("none"))     #Options: each - every plot will be printed in one pdf(Nodes) and one png (edges)
                                                           #one - all plots will be in one pdf(Nodes) and one pmg(edges)
                                                          #none - nothing will be saved
dev.off()

toc()
##cluster algos
# cfg <- cluster_label_prop(g)
# cfg <- cluster_fast_greedy(g)
# cfg <- cluster_louvain(g)
# cfg<-cluster_infomap(g)
# cfg<-cluster_walktrap(g)
# cfg<-cluster_spinglass(g)
# cfg <- cluster_edge_betweenness(g) #ACHTUNG DAUERT EWIG

dev.off()

tic("Calculate Heatmap")
cluster_informationen<-heatmap_clustered(igraph_object=igraph_object,
                                         cluster_algo=c("auto"),                                        #if set to auto optimal cluster algo will be chosen
                                         layout_for_network=layout,
                                         iterations=TRUE,                                               #iterations will produce a stable result but will take a little longer
                                         no_of_iterations=100,                                          #number of iterations cluster_louvain
                                         max_cluster_count_per_gene=10,                                  #no of clusters which a gene is allowd to be advised to before putting it into waste cluster
                                         min_cluster_size=15,                                            # min size of cluster to be shown in heatmap
                                         desicion_to_see_plot=F,                                     #whether to see the network plotted as well
                                         print_to_pdf=T,                                           #whether or not to save heatmap
                                         name_for_pdf=c("Heatmap_greater_clusters"),
                                         average = "mean")

summary[["heatmap_clustered"]]<-cluster_informationen$summary
clustered_heatmap_data<-cluster_informationen$heatmap_df

cluster_information<-cluster_informationen$color_cluster_min_size
View(cluster_information)

dev.off()
toc()

tic("Cluster_profiler")
############################################################################################################
##Cluster Profiler
##########################################################################################################
mart<-useMart("ensembl")
mart<-listDatasets(mart) # in dataframe names of available datasets attention only functions with human / mouse
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # load both independently from your organism
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl") # load both independently from your organism

clusterprofiler_results<-clusterprofiler_autoCena(cluster_to_check=c("greenyellow"),    #cluster to check can be either a certain cluster or all clusters
                                                  group=c("Disease"))         #condition you wanna check must be a column name of your info_Dataset

summary[["clusterprofiler_results"]]<-clusterprofiler_results$summary
clusterprofiler_results$summary<-NULL
toc()

