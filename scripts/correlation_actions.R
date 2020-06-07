#matrix operation valid for small dataset
#the result is a list of 3 matrices, (r) coefficient, (n) number of samples and (P) values
#it is not a data frame but a list of 3 df/matrices

correlation_matrix=Hmisc::rcorr(as.matrix(dd2))

#calculate number of combinations possible in pairwise gene vs gene comparisons (n*n-1/2)
#adj_p_value=FALSE
#if(adj_p_value==TRUE) {combinations = (nrow(correlation_matrix$P)*(nrow(correlation_matrix$P)-1))/2} else {combinations=1}

#do the bonferroni correction for the matrix p vals (adj pval = pval *number of comparisons)
#correlation_matrix$P= correlation_matrix$P*combinations
#create a 4 column data frame node1(gene1),node2(gene2), rval, pval(adj)

correlation_df = as.data.frame(t(combinat::combn(row.names(correlation_matrix$r), m=2)))

# write.table(correlation_df,
#             paste0(mainDir,"/reference_files/edge_list.csv"),
#             row.names = F,
#             col.names = T,
#             sep=",")
#correlation_df= read.csv(paste0(mainDir,"/reference_files/edge_list.csv"), header=T)

correlation_df$rval = correlation_matrix$r[lower.tri(correlation_matrix$r)]
correlation_df$pval = correlation_matrix$P[lower.tri(correlation_matrix$P)]

#retain rows which have pval (adj) < 0.05, and correlations above 0
correlation_df_filt = correlation_df[correlation_df$pval<0.05 & correlation_df$rval>0,]


#range of cutoff min to max (correlation)
range_cutoff<-seq(from = min_corr , to = max(correlation_df$rval) , length.out = range_cutoff_length)
range_cutoff<-round(range_cutoff, 3)
if(length(range_cutoff)>range_cutoff_length) {
  range_cutoff=range_cutoff[1:range_cutoff_length]
} else {
  range_cutoff=range_cutoff
}