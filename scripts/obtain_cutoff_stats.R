#@shobhit inform user about the number of cutoff(s) tested, the actual length (might be less than specified owing to rounding off to 3 digits)
#@following function should be parralelized

cutoff_prep=function(cutoff, corrdf_r, print.all.plots){
  
  # filtmatrix =  ifelse(correlation_matrix$r >cutoff, 1,0)
  # filtgraph = igraph::graph_from_adjacency_matrix(filtmatrix)
  
  ##define function for calculating graph statistics
  ##using only slightly modified version of the original algo
  #@shobhit edit function can do it on the filtered matrix
  rsquaredfun = function(graph_df, cutoff){
    #igraph object
    filt_igraph = igraph::graph_from_data_frame(graph_df,directed = F, vertices = NULL)
    #number of nodes
    num_nodes=igraph::vcount(filt_igraph)
    #number of edges
    num_edges=igraph::gsize(filt_igraph)
    #number of components
    graph_components <- igraph::components(filt_igraph)
    num_networks = graph_components$csize[graph_components$csize >= min_nodes_number_for_network] %>% 
      length()
    
    ##calculate stats
    d = igraph::degree(filt_igraph, mode = "all")
    dd = degree.distribution(filt_igraph, mode = "all", cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    
    
    
    if(length(probability)==0){
      R.square<-0
    }else{
      forplot = data.frame(probability=probability, degree=degree)
      reg = lm(log(probability) ~ log(degree))
      cozf = coef(reg)
      power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
      alpha = -cozf[[2]]
      R.square = summary(reg)$r.squared
      
     # print(paste("Alpha =", round(alpha, 3)))
     # print(paste("R square =", round(R.square, 3)))
     # plot
      
      #pdf("Degree_distribution_plot.pdf")
       # dd_plot= plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       #       col = 1, main = "Degree Distribution")
       # dd_plot + curve(power.law.fit, col = "red", add = T, n = length(d))
      #dev.off()
      
      if(print.all.plots) {
        
        degree_distribution_wd<-paste0("dir_DegreeDistribution_", topvar_genes)

        if(!degree_distribution_wd %in% list.dirs(working_directory)) {
          dir.create(paste0(working_directory,degree_distribution_wd))}

        dd_plot = ggplot(forplot,aes(x=log(degree), y= log(probability))) +
          geom_point() +
          geom_smooth(method="lm") +
          theme_bw()
        
        ggsave(path = paste0(working_directory,degree_distribution_wd),
               filename = paste0("Degree_distribution_plot_", cutoff, ".pdf"),
               dd_plot, device = cairo_pdf)
        
        
      }
      
       
    }
    
    
    output =  data.frame(R.squared=R.square,
                         degree=degree,
                         Probs=probability,
                         cutoff=cutoff,
                         no_edges=num_edges,
                         no_nodes=num_nodes,
                         no_of_networks=num_networks)
    
    return(output)
    
  }
  
  
  ###filter correlations above the cutoff
  filteredmatrix = corrdf_r[corrdf_r$rval>cutoff,]
  
  ##create expected df with initial values
  output=data.frame(R.squared=0,
                    degree=0,
                    Probs=0,
                    cutoff=cutoff,
                    no_edges=0,
                    no_nodes=0,
                    no_of_networks=0)
  
  #create a switch, if number of rows =0 then just output the initial output else do calculations
  #@shobhit check if the length of the nrow zero df and gtzero df are required to be of the same length
  rownums = ifelse(nrow(filteredmatrix)==0, "zero_rows", "gtzero")
  
  switch(rownums, zero_rows ={output},
         gtzero={rsquaredfun(graph_df=filteredmatrix, cutoff = cutoff)})
  
}


##what value of cutoff gives me ...


