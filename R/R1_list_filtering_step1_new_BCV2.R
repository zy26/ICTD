
BCV_ttest2_NA<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
{
  print("call new BCV test2 function !!!")
  x<-data0
  fff_cc<-c()
  for(kk in 1:rounds)
  {
    cv_result <- cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
    fff_cc<-rbind(fff_cc,cv_result$msep)
  }
  pp<-c()
  ddd<-apply(fff_cc,2,mean,na.rm=T)
  #print("bcv ttest2 ddd:")
  print(ddd)
  
  ddd<-ddd/sum(ddd)
  for(kk in 1:(ncol(fff_cc)-1))
  {
    pp_c<-1
    if(mean(fff_cc[,kk],na.rm=T)>mean(fff_cc[,kk+1],na.rm=T))
    {
      if(ddd[kk]>msep_cut)
      {
        pp_c<-t.test(fff_cc[,kk],fff_cc[,kk+1])$p.value
      }
    }
    pp<-c(pp,pp_c)
  }
  #boxplot(fff_cc)
  return(pp)       
}



R1_list_filtering_step1_new_BCV2 <- function (list_c2, data_CORS_cancer, max_cut = 20, cutn0 = 20, 
          cut10 = 0.8, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GPL570) 
{
  tg_1_rank_markers0 <- list_c2[[1]]
  tg_m_names0 <- list_c2[[2]]
  RMSE_on_CORS_cancer_c <- RMSE_row(data_CORS_cancer)
  tg_marker_lists <- list()
  for (i in 1:length(list_c2[[1]])) {
    ccc <- c(names(list_c2[[1]])[i], rownames(list_c2[[1]][[i]]))
    tg_marker_lists[[i]] <- ccc
  }
  names(tg_marker_lists) <- names(list_c2[[1]])
  pp_all <- c()
  for (i in 1:length(tg_marker_lists)) {
    pp <- sum(BCV_ttest2_NA(data_CORS_cancer[tg_marker_lists[[i]], 
                                          ], maxrank0 = 20, msep_cut = 0.01) < 0.001)
    pp_all <- c(pp_all, pp)
  }
  pp_R1_marker_list_f1 <- clean_rank1_module_new(data_CORS_cancer, 
                                                 tg_marker_lists, pp_all, st0 = 6)
  pp_R1_marker_list_f1.5 <- cut_modules(pp_R1_marker_list_f1, 
                                        cutn = cutn0)
  stat_p1 <- list()
  for (i in 1:length(pp_R1_marker_list_f1.5)) {
    tg_gene_c <- pp_R1_marker_list_f1.5[[i]][-1]
    stat_p1[[i]] <- tg_1_rank_markers0[[names(pp_R1_marker_list_f1.5)[i]]][tg_gene_c, 
                                                                           ]
  }
  names(tg_marker_lists) <- names(pp_R1_marker_list_f1.5)
  tg_genes_all <- c()
  for (i in 1:length(pp_R1_marker_list_f1.5)) {
    tg_genes_all <- c(tg_genes_all, pp_R1_marker_list_f1.5[[i]])
  }
  tg_genes_all <- unique(sort(tg_genes_all))
  print("filter 1 and stat done!")
  R1_markers_f1 <- pp_R1_marker_list_f1.5
  cut1 <- cut10
  ccc <- compute_min_jaccard(R1_markers_f1)
  ccc0 <- ccc > cut1
  stat_cc <- c(1:nrow(ccc0))
  names(stat_cc) <- 1:nrow(ccc0)
  for (i in 1:nrow(ccc0)) {
    for (j in 1:ncol(ccc0)) {
      if ((i < j) & (ccc0[i, j] > 0)) {
        nn <- max(i, j)
        stat_cc[which(stat_cc == i)] <- nn
        stat_cc[which(stat_cc == j)] <- nn
      }
    }
  }
  table(stat_cc)
  tg_ccc <- unique(stat_cc)
  R1_marker_list_f2 <- list()
  N <- 0
  for (i in 1:length(tg_ccc)) {
    N <- N + 1
    tg_ids <- as.numeric(names(stat_cc)[which(stat_cc == 
                                                tg_ccc[i])])
    ccc <- c()
    for (j in 1:length(tg_ids)) {
      ccc <- c(ccc, R1_markers_f1[[tg_ids[j]]])
    }
    ccc <- unique(ccc)
    R1_marker_list_f2[[N]] <- ccc
  }
  R1_marker_list_f2.5_stat <- rank_based_module_sorting(data_CORS_cancer, 
                                                        R1_marker_list_f2, IM_id_list, immune_cell_uni_table = immune_cell_uni_table)
  R1_marker_list_f2.5 <- R1_marker_list_f2.5_stat[[1]]
  pp_all <- c()
  for (i in 1:length(R1_marker_list_f2.5)) {
    pp <- sum(BCV_ttest2_NA(data_CORS_cancer[R1_marker_list_f2.5[[i]], 
                                          ], maxrank0 = 20, msep_cut = 0.01) < 0.001)
    pp_all <- c(pp_all, pp)
  }
  pp_R1_marker_list_f3 <- clean_rank1_module(data_CORS_cancer, 
                                             R1_marker_list_f2.5, pp_all, st0 = 6)
  pp_R1_marker_list_f3.5 <- cut_modules(pp_R1_marker_list_f3, 
                                        cutn = cutn0)
  R1_marker_list_f3.5_stat <- rank_based_module_sorting(data_CORS_cancer, 
                                                        pp_R1_marker_list_f3.5, IM_id_list, immune_cell_uni_table = immune_cell_uni_table)
  print("filter 2 done!")
  ccc <- c()
  nn <- c()
  for (i in 1:length(pp_R1_marker_list_f3.5)) {
    ccc0 <- c()
    for (j in 1:length(IM_id_list)) {
      if (length(IM_id_list[[j]]) > 1) {
        cc0 <- apply(immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]], 
                                           IM_id_list[[j]]], 1, sum)/sum((1/(1:length(IM_id_list[[j]]))))
      }
      else {
        cc0 <- immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]], 
                                     IM_id_list[[j]]]
      }
      ccc0 <- cbind(ccc0, cc0)
    }
    colnames(ccc0) <- names(IM_id_list)
    ddd <- apply(ccc0, 2, mean)
    ccc <- rbind(ccc, ddd)
    nn <- c(nn, colnames(ccc0)[which(ddd == max(ddd))[1]])
  }
  rownames(ccc) <- nn
  cell_enrich_stat <- ccc
  names(pp_R1_marker_list_f3.5) <- nn
  rrr <- rep(1, length(pp_R1_marker_list_f3.5))
  Filter_1_result_list <- list(pp_R1_marker_list_f1, R1_marker_list_f2, 
                               R1_marker_list_f3.5_stat, pp_R1_marker_list_f3.5, rrr, 
                               cell_enrich_stat)
  names(Filter_1_result_list) <- c("R1_marker_list_f1", 
                                   "R1_marker_list_f2", "R1_marker_list_f3.5_stat", 
                                   "R1_marker_list_f3.5", "R1_marker_list_rank", 
                                   "R1_marker_list_f3.5_cell_enrich_stat")
  return(Filter_1_result_list)
}