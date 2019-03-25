#######################################CoCena²#########################
### R version check
if(grepl(c("3.4.4"),R.version.string)){
  print("Welcome to CoCena²")
}else{
  print("You can try but there can be troubles with other versions")
}
######################################################################

#####Set working directory
setwd("E:/RNA-Seq/Wataru_COPD/Data/output/kallisto/")
mainDir <- getwd()
subDir <- c("CoCena")
dir.create(file.path(mainDir , subDir))
setwd(file.path(mainDir , subDir))
originalwd <- getwd()

# load packages

packages_to_load <- c("igraph" , "plotly" ,"ggplot2" ,"bnstruct", "gridExtra" , "Hmisc" , 
                      "RColorBrewer" , "gtools" , "rlist" , "reshape2" , "gplots" , "moduleColor" , 
                      "NMF" , "rlist" , "ggpubr","RJSONIO","httr","stringr","XML","DOSE","clusterProfiler",
                      "devtools","clues","bnstruct","fields","curl","httpuv","intergraph","bnlearn","fields","viridis","foreign")

lapply(packages_to_load , require , character.only=TRUE)
install_github('cytoscape/cytoscape-automation/for-scripters/R/r2cytoscape')

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
library(readr)
library(circlize)

library(foreign)
library(httpuv)

#install.packages("devtools")
library(devtools)
#install_github("kassambara/r2excel")
library(r2excel)
library(grid)
library(clusterProfiler)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ReactomePA")
library(ReactomePA)
#source("https://bioconductor.org/biocLite.R")
#biocLite("pcaGoPromoter")
library(pcaGoPromoter)
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("biomaRt")

library(biomaRt)
library(reshape2)
library(dplyr)
require("AnnotationDbi")
library(org.Hs.eg.db) 
library(org.Mm.eg.db) 
library(dendextend) 

source("E:/owncloud/backup/R/Cocena_functions.R")


#####load in your data

###############################################################################################
#make sure that in your dataset, 
#gene names are as rownames and samples are as column names!
###############################################################################################
removedbatch_dds_vst <- limmaBatchEffectRemoval(input=dds_vst,
                                                modelfactor = "condition",
                                                batchfactor = c("sex"),
                                                batchfactor_2 = "date.of.library")


Dataset_1 <- removedbatch_anno_log2
Dataset_1 <- Dataset_1[,c(1:(ncol(Dataset_1)-3))]
Dataset_1 <- Dataset_1[,c(-(ncol(Dataset_1)-1))]
Dataset_1_symbol <- Dataset_1[!duplicated(Dataset_1[,ncol(Dataset_1)]), ]
Dataset_1_symbol <- Dataset_1_symbol[complete.cases(Dataset_1_symbol), ]

row.names(Dataset_1_symbol) <- Dataset_1_symbol[,ncol(Dataset_1)]
Dataset_1 <- Dataset_1_symbol[,c(-(ncol(Dataset_1)))]
#Dataset_1 <- Dataset_1[,c(4, 10, 16, 2, 8, 14, 6, 12, 18, 3, 9, 15, 1, 7, 13, 5, 11, 17)]

View(Dataset_1)
#Dataset_1<-t(Dataset_1)
###############################################################################################
#make sure that sample names are in rownames
#MAKE sure that there is no " - "in merged names - cytoscape cannot handle these
#Make sure that there is no space in your names- cytoscape cannot handle these
######################################################################################
info_Dataset <- sample_table

#info_Dataset <- info_Dataset[order(info_Dataset$merged),] 
rownames(info_Dataset)<-info_Dataset$ID
info_Dataset$ID<-NULL

info_Dataset <-info_Dataset[,c(1:2)]
colnames(info_Dataset) <- c("merged","conditions")


#info_Dataset <- info_Dataset[!rownames(info_Dataset) %in% c("6517"), ]


#info_Dataset <- as.data.frame(info_Dataset[c(4, 10, 16, 2, 8, 14, 6, 12, 18, 3, 9, 15, 1, 7, 13, 5, 11, 17),])

View(info_Dataset)

#####################################################################################
#set the organism your data is from
organism=c("human") #//c("Mouse")
####################################################################################
##load additional needed Information!

library(clusterProfiler)
TF_list <- read.delim("E:/RNA-Seq/Viemann_lnRNAs/lncRNA_Bachelor_Marie/Data/TFcat.txt" , header = TRUE , check.names = FALSE)
gmtfile <- clusterProfiler::read.gmt("E:/Paper/DeNardo-WGCNA/h.all.v5.2.symbols.gmt")
epigenetic_modulators <- read_delim("E:/Listen/epigenetic_modulators.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
GO_immune_response <- xlsx.readFile("E:/RNA-Seq/Viemann_lnRNAs/lncRNA_Bachelor_Marie/Data/GO_immune_response.xlsx", sheetIndex =1)
GO_metabolic_process <- xlsx.readFile("E:/RNA-Seq/Viemann_lnRNAs/lncRNA_Bachelor_Marie/Data/GO_metabolic_process.xlsx", sheetIndex =1)
GO_signalling <- xlsx.readFile("E:/RNA-Seq/Viemann_lnRNAs/lncRNA_Bachelor_Marie/Data/GO_signalling.xlsx", sheetIndex =1)
####################################################################################


#Save original data to new variable which could be changed during process
#original Dataset can always be found under Dataset_1
original_data <- Dataset_1


#####################################################################################
##filter data if wanted
filter=T
#####################################################################################
#1) on transcription factors
original_data <- original_data[rownames(original_data) %in% TF_list$Mouse , ] ## HUman or mouse

#2) any other gene_list
gene_list <- read.table("E:/RNA-Seq/Creld/CoCena_not_switched_woProt/union_only_significant_sinlge_comp_0.01_Anna.txt" , header = TRUE , check.names = FALSE)

original_data <- original_data[rownames(original_data) %in% gene_list$ID , ] ## HUman or mouse

original_data <- Dataset_1


#3) topX
dd1 = original_data[order(apply(original_data, 1, var),decreasing = T),]
original_data <- head(dd1,2000)




original_data_search <- original_data
original_data_search$symbol <- row.names(original_data)

####################################################################################
####################################################################################
collectGarbage()

#####
#first visualisation of your data
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

#original_data <- original_data[apply(original_data,1,function(x) !sd(x)==0),]

#original_data <- original_data[1:4000,]

pheatmap::pheatmap(original_data, 
                   cluster_row = T,
                   cluster_cols = T,
                   color = col.pal,
                   scale = c("row"),
                   show_rownames = F,
                   annotation_col = info_Dataset,
                   main = c("Heatmap: Gene Expression Data - whole Dataset "))

####################################################################################
#There will be a function collecting all set variables values in this list "summary"
###
summary<-list()
####################################################################################
####################################################################################


###Start CoCena²
correlation_df<-correlation(type_of_correlation="pearson",            #pearson / spearman
                            pVal_treshold=0.05,                       #pValue Treshold for calculated correlation
                            save_raw_data=FALSE)                      #save correlation data in .txt

summary[["correlation"]]<-correlation_df$summary
correlation_df<-correlation_df$raw_data


####deciding on cutoff
library(igraph)
cutoff_options<-cutoff_visualisation(correlation_df=correlation_df,
                                     min_corr=0.7,                    # min correlation 
                                     range_cutoff_length=60,          #No.of cutoffs tested (all higher than min_cor)
                                    min_no_for_cluster=10)            #how many nodes are minimum for a network 

summary[["cutoff_visualisation"]]<-cutoff_options$summary

####################################################################################
##Look at plot
#cutoff_options$Plotly_object
#slide through cutoffs and decide on one
#noch händisch einzutragen
####################################################################################
cutoff_options$Cutoff_df
chosen_cutoff<-0.750 #top1000 17 samples
chosen_cutoff<-0.814 #top5000 17 samples
chosen_cutoff<-0.883 #top2000 15 samples
chosen_cutoff<-0.882 #top2000 15 samples sex + lib
chosen_cutoff<-0.87 #top5000 15 samples
chosen_cutoff<-0.888 #top10000 15 samples


###Change to you chosen cutoff!!
####################################################################################
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
                            layout=c("layout_with_fr"),                          #choose anything       if you tested layouts change to layout_on_grid to fasten things up
                            min_nodes_number_for_network=20,
                            show_HC= TRUE)                                       #see HC after cutting whole Dataset
summary[["plot_network"]]<-return_list$summary
igraph_object<-return_list[["graph_object"]]

######Scale free Topology check
fit_power_law(igraph_object)

####layout can be chosen from
layout<-return_list[["layout"]]

####################################################################################
######if test_layout was run
names(layouts_on_list_coord)# names of your layout options
layout<-layouts_on_list_coord$layout_with_fr
####################################################################################
#GFC calculation
GFC_all_genes<-GFC_calculation(normdata=original_data,
                               group=c("merged"),       #column name from info_Dataset
                               data_in_log=T,          #is the inloaded data at some point put in log (e.g. DeSeq's rlog)
                               range_GFC=2)              #if calculated GFC value is over range_GFC it is 
                                                           #set to range_GFC value due to nice visulaisation
summary[["GFC_calculation"]]<-GFC_all_genes$summary
GFC_all_genes<-GFC_all_genes$GFC_all_genes
#####error !!!!!
plot_GFC_networks(igraph_object=igraph_object,
                  print_to_pdf=F,                      #GFC plots will be plotted in pdf - ATTENTION when further work will be done in Corel Draw not recommended
                  print_edges_png_nodes_pdf=c("none"))     #Options: each - every plot will be printed in one pdf(Nodes) and one png (edges)
                                                           #one - all plots will be in one pdf(Nodes) and one pmg(edges)
                                                          #none - nothing will be saved
dev.off()


##cluster algos
# cfg <- cluster_label_prop(g)
# cfg <- cluster_fast_greedy(g)
# cfg <- cluster_louvain(g)
# cfg<-cluster_infomap(g)
# cfg<-cluster_walktrap(g)
# cfg<-cluster_spinglass(g)
# cfg <- cluster_edge_betweenness(g) #ACHTUNG DAUERT EWIG
##

#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)
cluster_informationen<-heatmap_clustered(igraph_object=igraph_object,
                                       cluster_algo=c("cluster_louvain"),                                        #if set to auto optimal cluster algo will be chosen
                                       layout_for_network=layout,
                                       iterations=TRUE,                                               #iterations will produce a stable result but will take a little longer
                                       no_of_iterations=100,                                          #number of iterations
                                       max_cluster_count_per_gene=10,                                  #no of clusters which a gene is allowd to be advised to before putting it into waste cluster
                                       min_cluster_size=5,                                            # min size of cluster to be shown in heatmap
                                       desicion_to_see_plot=TRUE,                                     #whether to see the network plotted as well
                                       desicion_to_save_plot=FALSE,                                   #whether or not to save network 
                                       name_for_pdf_plot=c("Network_greater_clusters_clustered"),
                                       print_to_pdf=T,                                           #whether or not to save heatmap
                                       name_for_pdf=c("Heatmap_greater_clusters_clustered"),
                                       average = "mean")


summary[["heatmap_clustered"]]<-cluster_informationen$summary
clustered_heatmap_data<-cluster_informationen$heatmap_df

cluster_information<-cluster_informationen$color_cluster_min_size



cluster_to_excel(sample_cluster = T)


########
##unzufrieden mit der Heatmap?
##sie haben zu viele cluster die den selben trend über alle conditionen aufzeigen?
##sie möchten selber bestimmen wie viele cluster sie haben wollen ?
##dann benutzen sie folgende funktion
merged_cluster_data<-merge_cluster(clustered_heatmap_data,cluster_wanted = 8)

dev.off()

#####



plot(hclust(dist(clustered_heatmap_data)), hang = -1, cex = 0.6)

height=1.4 #top5000
height=0.9 #top2000
height=1.0 #top2000
height=1.1 #top2000 sex + lib


abline(h=height,col="red")

merge_cluster<- merge_cluster_Thomas(height=height,cluster_cols=T)
  

summary[["merge_cluster"]]<-merge_cluster$summary
cluster_information<-merge_cluster$cluster_information_new


cluster_to_excel(sample_cluster = T)


###################################################################################
#Circos Plot for all clusters

y_compare <- compareClusterGO() # This takes very long. This step can be skipped,
                                # then you just need to set the y_comp parameter
                                # in clusterCircos() to F.

heatmap_vec <- c("GFC_COPD_2",  "GFC_Control" ,"GFC_COPD_3"  ,"GFC_COPD_4")

dev.off()

clusterCircos(hm_vec = TRUE, seg = 100, range = c(-1, 0, 1), y_comp = T)




################################to cytoscape########################################
###CYTOSCAPE MUST BE OPEN
toCytoscape_all(cluster_information=cluster_information)
##getting the layout coordinates
layout<-as.matrix(LayoutfromCytoscape()) ## rerun from pltGFCnetworks

####################################################################################

#testing out different cutoffs
search_for_good_cutoff(data=correlation_df,
                       min_nodes_number_for_network=10,            #min no of nodes connected to count as network
                       show_network=TRUE,                          #show plotted network
                       layout=layout_with_kk,                      #choose layout algorithm
                       cluster_algo= cluster_louvain,              # choose cluster algorithm
                       max_cluster_count_per_gene=8,               #no of clusters which a gene is allowd to be advised to before putting it into waste cluster
                       min_cluster_size=10,                         #min size of cluster to be shown in heatmap
                       abberation=0.1,                             #abberation from before chosen cutoff in percent , can be set to 0
                       no_cutoffs_tested=4,                        # no of cutoffs you want to test
                       cutoffs_to_test=c("0.5","0.8","0.2"))      #set abberation to 0 if you want certain cutoffs tested
####output ?! maybe list with the igraph objects ?!


unique(cluster_information$try)


############################################################################################################
##Cluster Profiler
##########################################################################################################
mart<-useMart("ensembl")
mart<-listDatasets(mart) # in dataframe names of available datasets attention only functions with human / mouse
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") # load both independently from your organism
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl") # load both independently from your organism

clusterprofiler_results<-clusterprofiler_autoCena(cluster_to_check=c("all"),    #cluster to check can be either a certain cluster or all clusters
                                                  group=c("merged"))         #condition you wanna check must be a column name of your info_Dataset

summary[["clusterprofiler_results"]]<-clusterprofiler_results$summary
clusterprofiler_results$summary<-NULL
###
#following creates Dataframe where each gene gets its annotation and a color depending in what analysis it was called
# color_cluster_info<-cluster_information_to_color(vertex_attributes=cluster_information,        # Dataframe with every gene and its belonging cluster
#                                                  clusterprofiler_results=clusterprofiler_results) # cluster profiler results to add
# summary[["cluster_information_to_color"]]<-color_cluster_info$summary
# color_cluster_info<-color_cluster_info$color_saved
# color_cluster_info<-as.data.frame(color_cluster_info)

igraph_object_all<-igraph_object   # big network is saved in  igraph_object_all // parent network



####################################################
#Plot Gene on network###############################
####################################################
####################################################
library(data.table)
kegg_terms <- file("E:/RNA-Seq/Wataru_COPD/Data/output/kallisto/c2.all.v6.2.symbols.gmt")
kegg_terms <- read.gmt(kegg_terms)
kegg_terms <- kegg_terms[c("gene", "ont")]
kegg_terms <- data.frame(lapply(kegg_terms, as.character), stringsAsFactors=FALSE)

kegg_terms_list <- c("LIPID","MIGRATION","ENDOCYTOSIS","GLYCOLYSIS")
term = "LIPID"
term = "MIGRATION"
term = "ENDOCYTOSIS"
term = "GLYCOLYSIS"


for(geneset_name in kegg_terms_list){

  goi <- kegg_terms[kegg_terms$ont %like% geneset_name, ]
  plot_genes_on_networks(igraph_object=igraph_object,geneset_name = geneset_name, geneset = goi$gene,print_to_pdf=T)
}




#Metabolism

reactome_terms <- file("E:/Listen/GeneSets/reactome.v5.2.symbols_mouse.gmt")
reactome_terms <- read.gmt(reactome_terms)
reactome_terms <- reactome_terms[c("gene", "ont")]
reactome_terms <- data.frame(lapply(reactome_terms, as.character), stringsAsFactors=FALSE)

reac_goi_1 <- intersect(reactome_terms[reactome_terms$ont =="REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS",]$gene,cluster_information$Gene)
reac_goi_2 <- intersect(reactome_terms[reactome_terms$ont =="REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",]$gene,cluster_information$Gene)
reactome <- union(reac_goi_1,reac_goi_2)


hallmark_terms <- file("E:/Listen/GeneSets/msigdb.v5.2.symbols_mouse.gmt")
hallmark_terms <- read.gmt(hallmark_terms)
hallmark_terms <- hallmark_terms[c("gene", "ont")]
hallmark_terms <- data.frame(lapply(hallmark_terms, as.character), stringsAsFactors=FALSE)

hall_goi_1 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_GLYCOLYSIS",]$gene,cluster_information$Gene)
hall_goi_2 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_OXIDATIVE_PHOSPHORYLATION",]$gene,cluster_information$Gene)
hall_goi_3 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_FATTY_ACID_METABOLISM",]$gene,cluster_information$Gene)

hallmark <-union(union(hall_goi_1,hall_goi_2),hall_goi_3)


kegg_terms <- file("E:/Listen/GeneSets/kegg.v5.2.symbols_mouse.gmt")
kegg_terms <- read.gmt(kegg_terms)
kegg_terms <- kegg_terms[c("gene", "ont")]
kegg_terms <- data.frame(lapply(kegg_terms, as.character), stringsAsFactors=FALSE)

kegg_goi_1 <- intersect(kegg_terms[kegg_terms$ont =="KEGG_CITRATE_CYCLE_TCA_CYCLE",]$gene,cluster_information$Gene)
kegg_goi_2 <- intersect(kegg_terms[kegg_terms$ont =="KEGG_OXIDATIVE_PHOSPHORYLATION",]$gene,cluster_information$Gene)
kegg <- union(kegg_goi_1,kegg_goi_2)


union_metabolism <- union(union(kegg,hallmark),reactome)

plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "Metabolismn_kegg", geneset = kegg, print_to_pdf = T)
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "Metabolismn_hallmark", geneset = hallmark, print_to_pdf = T)
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "Metabolismn_reactome", geneset = reactome, print_to_pdf = T)



#MTOR signaling

hallmark_terms <- file("E:/Listen/GeneSets/msigdb.v5.2.symbols_mouse.gmt")
hallmark_terms <- read.gmt(hallmark_terms)
hallmark_terms <- hallmark_terms[c("gene", "ont")]
hallmark_terms <- data.frame(lapply(hallmark_terms, as.character), stringsAsFactors=FALSE)

hall_goi_1 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_PI3K_AKT_MTOR_SIGNALING",]$gene,cluster_information$Gene)
hall_goi_2 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_MTORC1_SIGNALING",]$gene,cluster_information$Gene)

hallmark <-union(hall_goi_1,hall_goi_2)
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "MTOR signaling_hallmark", geneset = hallmark, print_to_pdf = T)



#Cell cycle
hall_goi_1 <- intersect(hallmark_terms[hallmark_terms$ont =="KEGG_DNA_REPLICATION",]$gene,cluster_information$Gene)
hall_goi_2 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_CELL_CYCLE_CHECKPOINTS",]$gene,cluster_information$Gene)
hall_goi_3 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_MITOTIC_SPINDLE",]$gene,cluster_information$Gene)
hall_goi_4 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_G1_PHASE",]$gene,cluster_information$Gene)
hall_goi_5 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_G1_S_TRANSITION",]$gene,cluster_information$Gene)
hall_goi_6 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_P53_PATHWAY",]$gene,cluster_information$Gene)
hall_goi_7 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_MITOTIC_G2_G2_M_PHASES",]$gene,cluster_information$Gene)
hall_goi_8 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_NUCLEOTIDE_EXCISION_REPAIR",]$gene,cluster_information$Gene)

cell_cycle <- Reduce(union, list(hall_goi_1,hall_goi_2,hall_goi_3,hall_goi_4,hall_goi_5,hall_goi_6,hall_goi_7,hall_goi_8))
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "Cell cycle", geneset = cell_cycle, print_to_pdf = T)







#Wnt/b-catenin signaling

hall_goi_1 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_WNT_BETA_CATENIN_SIGNALING",]$gene,cluster_information$Gene)
hall_goi_2 <- intersect(hallmark_terms[hallmark_terms$ont =="REACTOME_SIGNALING_BY_WNT",]$gene,cluster_information$Gene)

wnt <- Reduce(union, list(hall_goi_1,hall_goi_2))
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "WNT", geneset = wnt, print_to_pdf = T)



#Inflammatory response

hall_goi_1 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_INFLAMMATORY_RESPONSE",]$gene,cluster_information$Gene)
hall_goi_2 <- intersect(hallmark_terms[hallmark_terms$ont =="HALLMARK_TNFA_SIGNALING_VIA_NFKB",]$gene,cluster_information$Gene)
hall_goi_3 <- intersect(hallmark_terms[hallmark_terms$ont =="KEGG_JAK_STAT_SIGNALING_PATHWAY",]$gene,cluster_information$Gene)

wnt <- Reduce(union, list(hall_goi_1,hall_goi_2,hall_goi_3))
plot_genes_on_networks(igraph_object=igraph_object,geneset_name = "Inflammatory response", geneset = wnt, print_to_pdf = T)


####################################################
####################################################

#############################################################################################################
#############################################################################################################
########################################## Intraclusteral Analysis ##########################################
############################################################################################################
#cluster names you can chose from:
unique(cluster_information$try)

list_of_results<-plot_single_cluster(igraph_object=igraph_object_all,               #parent network 
                                     cluster_name=c("turquoise"),        #cluster you chosen
                                     top_percentage_for_hubs=0.25,                  #percentage of how many genes is allowed to be hubs
                                     allowed_edges=20,                              #allowed edges if gene is counted as hub
                                     allowed_edges_between_hubs=2,                  #allowed edges between hubs
                                     string_needs_to_be_redone=TRUE,
                                     string_treshold=500,                           #if testing different parameters you not always need to redo string -
                                     color_STRING=c("#c95555"),                     #color to mark edges stated in STRING
                                     color_edges=c("grey"),                         # color proposed edges /edges not found on STRING
                                     no_strings_to_be_string_hub=1,                   #String edges can bring forward new HUBS, 
                                                                                    #determine how many you want to have allowed - 
                                                                                    #percentage is depending on degree number
                                     label_all_TF=TRUE,                            #TRUE/FALSE is you want all Transcription factors labeld - 
                                                                                    #if FALSE only labelled if these fall under categorie labelled
                                     color_label_if_TF=c("#D3436EFF"),              # color of TF
                                     color_label_normal=c("#231151FF"),             # color of not TF
                                     width_label_edge=2,                            #width label edge to be more striking
                                     width_normal=1,                                # width normal edge
                                     percentage_named=0.7,                          # how many genes you want to be labeld depending an all paramteres before
                                     size_label_boxes=30)                           # adjust label boxes size normally between 20 and 30
summary[["I-GIN"]]<-list_of_results$summary
to_save<-list_of_results$summary[2:length(list_of_results$summary)]
to_save<-list.rbind(to_save)
summary[["Output_IGIN"]]<-list_of_results[["igraph_object_small"]]
igraph_object_small<-list_of_results[["igraph_object_small"]]
vertex_attributes_label<-list_of_results[["vertex_attributes_cluster"]]
edges_selected_label<-list_of_results[["edges_attributes"]]
############################################################################################################
#########################################Plot this beautiful thing##########################################
###########################################################################################################
 library(qgraph)

e <- igraph::get.edgelist(igraph_object_small,names = FALSE)

l <- qgraph::qgraph.layout.fruchtermanreingold(e,vcount=igraph::vcount(igraph_object_small),
                                               weights = igraph::edge.attributes(igraph_object_small)$weight,
                                       area=10*(igraph::vcount(igraph_object_small)^2*2),                          #play around with area and reoulse.rad do minimise overlap and
                                       repulse.rad=(igraph::vcount(igraph_object_small)^3.1),                        #achieve nice plotting
                                       niter = 1000)                                                             #the higher niter the longer it takes but also more trys to improve

#################
#IGIN
################
#pdf("Übelster_nicest_network.pdf", height = 10, width = 16)
#par(mar=c(5.1,4.1,4.1,2.1))

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


#############################search and label certain genes by hand########################################
add_gene_label(gene_to_add=c("Trp53"))



##############################################################################
###plot bayesian_based intresting genes with GFC indicating it's expression###
bayesian_based_subnetworks(igraph_object_small)

pieplot_coloring_plot(igraph_object=igraph_object_small, 
                      degree_to_color=FALSE,degree_to_color_percentage=0.1,
                      bayesian_to_color=TRUE,
                      resize_nodes=TRUE,
                      size_parent_pies=15,
                      size_children_pies=10,
                      size_others=1,
                      edge_width=2)

dev.off()
# bayesian_network_igraph<-graph_from_data_frame(bnnet2$arcs,directed = TRUE)
# bayesian_network_igraph
# V(bayesian_network_igraph)$name
# label_options_bayesian<-data.frame(Gene=V(bayesian_network_igraph)$name, label = " ", color="lightgreen",size=4,stringsAsFactors = FALSE)
# str(label_options_bayesian)
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"label"]<-c("Mmp14")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"label"]<-c("Runx3")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"label"]<-c("Nek6")
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"color"]<-"yellow"
# label_options_bayesian[label_options_bayesian$Gene %in% c("Mmp14"),"size"]<-15
# label_options_bayesian[label_options_bayesian$Gene %in% c("Runx3"),"size"]<-15
# label_options_bayesian[label_options_bayesian$Gene %in% c("Nek6"),"size"]<-15
# shortestPath<-shortest_paths(bayesian_network_igraph,from="Arhgap24",to="Runx3")
# 
# E(bayesian_network_igraph, path = shortestPath$vpath[[1]])$color<-"red"
# graph_from_edgelist()
# 
# plot(bayesian_network_igraph,
#      layout=layout_as_tree,
#      #vertex.color=label_options_bayesian$color,
#      # vertex.label=label_options_bayesian$label,
#      vertex.size=5,
#      edge.arrow.size=0.5)


######saving
mainDir <- cutoff_wd
subDir <- c("IGIN")
dir.create(file.path(mainDir , subDir))
setwd(file.path(mainDir , subDir))
IGIN_wd<-getwd()



write.table(to_save,paste0(IGIN_wd,"/",to_save["cluster_name",],".txt"),sep="\t")



#####LEGENDE eher semi

# #par(mar=c(5.1,0.5,4.1,0.5))
# colfunc <- colorRampPalette(c('#B1B1B2','#982D80FF'))
# graphics::legend(x=-2,y=1.5, legend = c("Edges","known interactions","unknown interactions"," ","Labels","transcription factor","no transcription factor"," ","Nodes"),
#        col = c("white","#c95555","grey","white","white","#D3436EFF","#231151FF", "white") , pch = c(15,15,15,15,15,15,15,15),
#        bty = "n", pt.cex = 2.2, cex = c(1,1,1,1,1,1,1), horiz = F , text.font = c(2,1,1,1,2,1,1,2))
# 
# 
# xl <- 1
# yb <- 0.4
# xr <- 1.5
# yt <- -0.1
# 
# graphics::rect(-1.93,-head(seq(yb,yt,(yt-yb)/20),-1),-1.88,-tail(seq(yb,yt,(yt-yb)/20),-1), col = colfunc(20), border = NA)
# 
# 
# 
# 
# mtext(c("receiving > regulating"," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," "," ","regulating > receiving"),side=4,at=-tail(seq(yb,yt,(yt-yb)/10),20),line=-62,las=2,cex=1)
# 
# 
# 
# dev.off()



#############################################################################################################
#################################################  summary  ##################################################
##############################################################################################################
#big list with included data
summary #can load some time !

#output most important set values during script
summarized_summary<-summary_wrapped_up(print_text = TRUE)

list.save(summary, "summary.rdata") # will save all data which is important to run the script, there is another script
# optimised for this summary data which can be passed along and reloaded on every other computer
#list.load("summary.rdata")


############################################################################################################
############################################################################################################
########################################### MARIE WORK #####################################################
############################################################################################################

BayesOutput <- BayesNet(clust.name =c("turquoise") ,
                        cutoff = chosen_cutoff , 
                        CPratio = 5 ,
                        PCratio = 3)

BayesOutput$BayesCompleteList$Mmp14.children
BayesOutput$BayesCompleteList$Mmp14.parents
BayesOutput$BayesCompleteList$Runx3.children
BayesOutput$BayesCompleteList$Runx3.parents
#################################################################
##testing genes

# test_genes<-c("Runx3")
# test_genes_expression<-Dataset_1[test_genes,]
# test_genes_expression<-rbind(test_genes_expression,as.character(info_Dataset$merged))
# test_genes_expression<-rbind(test_genes_expression,ifelse(grepl("TAM", info_Dataset$merged),c("TAM"),c("Control")))
# test_genes_expression<-as.data.frame(t(test_genes_expression))
# test_genes_expression$gene<-test_genes
# colnames(test_genes_expression)<-c("expression","Condition","group","Gene")
# # condition_tam<-data.frame(test_genes_expression[,grepl("TAM",info_Dataset$merged)],condition="condition",check.names=FALSE)
# #
# # control<-data.frame(test_genes_expression[,!(grepl("TAM",info_Dataset$merged))],condition="control",check.names=FALSE)
# test_genes_expression$expression<-as.numeric(as.character(test_genes_expression$expression))
# 
# test_genes_expression$group<-as.factor(test_genes_expression$group)
# test_genes_expression$Condition<-as.factor(test_genes_expression$Condition)
# Mmp14_test_genes_expression<-test_genes_expression
# test_genes<-c("Mmp14")
# test_genes_expression<-Dataset_1[test_genes,]
# test_genes_expression<-rbind(test_genes_expression,as.character(info_Dataset$merged))
# test_genes_expression<-rbind(test_genes_expression,ifelse(grepl("TAM", info_Dataset$merged),c("TAM"),c("Control")))
# test_genes_expression<-as.data.frame(t(test_genes_expression))
# test_genes_expression$gene<-test_genes
# colnames(test_genes_expression)<-c("expression","Condition","group","Gene")
# # condition_tam<-data.frame(test_genes_expression[,grepl("TAM",info_Dataset$merged)],condition="condition",check.names=FALSE)
# #
# # control<-data.frame(test_genes_expression[,!(grepl("TAM",info_Dataset$merged))],condition="control",check.names=FALSE)
# test_genes_expression$expression<-as.numeric(as.character(test_genes_expression$expression))
# 
# test_genes_expression$group<-as.factor(test_genes_expression$group)
# test_genes_expression$Condition<-as.factor(test_genes_expression$Condition)
# Runx3_test_genes_expression<-test_genes_expression
# 
# all<-rbind(Runx3_test_genes_expression,Mmp14_test_genes_expression)
# all$Gene<-as.factor(all$Gene)
# 
# all_without_PRE<-all[!(grepl("PRE",all$Condition)),]
# 
# ggplot(all,aes(x=group,y=expression,fill=Gene))+geom_boxplot()+geom_jitter()
# 
# dev.off()

#######

