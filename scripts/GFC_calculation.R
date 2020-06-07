GFC_calculation <- function(voi_id) {
  corresp_info = info_dataset[rownames(dd2)%in%rownames(info_dataset),]



  if(intersect(voi_id,colnames(info_dataset))%>% length()>0){

    print(paste("variables:",paste0(intersect(voi_id,colnames(info_dataset)), collapse = ","),
                "will be used as grouping variables,",
                "if these are different in # from voi_id then please check if",
                "the variables mentioned in voi_id are present in the metadata"))

    corresp_info$grpvar =purrr::pmap(corresp_info[intersect(voi_id,colnames(info_dataset))],
                                     paste, sep="-") %>% unlist()

  } else {

    print(paste("the first column in the metadata will be used as the grouping variable",
                "since the voi_id is not present in the metadata"))
    corresp_info$grpvar = info_dataset[,1]

  }


  norm_data_anno = merge(dd2, corresp_info["grpvar"], by="row.names",all.x=T)
  #colnames(norm_data_anno)[1] = "sample_names"
  #contains a column called Row.names (1st col)
  norm_data_anno = norm_data_anno[,-1]
  ##@shobhit original code must be parallelized by group
  #
  norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
  trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])
  #
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
  grplist=colnames(trans_norm)[-(ncol(trans_norm))]


  group_means= trans_norm[,"group_mean"]

  gfc_calc=function(grp){
    df1= trans_norm[,grp]
    df2=gtools::foldchange(df1, group_means) %>% ifelse(.>2.0, 2.0,.) %>% ifelse(.< (-2.0), -2.0,.) %>% as.data.frame()
    colnames(df2) = paste0("GFC_",grp)
    return(df2)
  }

  GFC_all_genes=do.call("cbind", lapply(grplist, gfc_calc))
  GFC_all_genes= round(GFC_all_genes,3)
  GFC_all_genes$Gene = rownames(GFC_all_genes)

  return(GFC_all_genes)


}
