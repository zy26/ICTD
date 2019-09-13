

print("test docker with Hello word!")

# Then, test the ICTD
library('ICTD')
#devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
library('EPIC')
library('readr')
library('nnls')

#data_bulk = GSE72056_diri_example[[1]]
#print(data_bulk[1:5,1:6])
# ictd_result <- ICTD(data_bulk)
# ictd_result[[2]]
#print(sessionInfo())


#----------------function part---------------
ictd_2_output <- function(ictd_prop, dataset.name)
{
  #2
  col_name <- matrix(colnames(ictd_prop),1,length(colnames(ictd_prop)))
  vv <- c()
  for(i in 1:length(col_name))
  {
    vv <- c(vv, rep(col_name[i], 8))
  }
  sample_col <- vv
  
  #3  
  celltype_col <- rep(c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts'), ncol(ictd_prop))
  #4
  prediction_col <- rep(runif(8,min=0.01,max=0.8), ncol(ictd_prop))  
  
  #bind
  output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)
  
  return(output_df)
}

ictd_2_output_real <- function(ictd_prop, dataset.name)
{
  #2
  col_name <- matrix(colnames(ictd_prop),1,length(colnames(ictd_prop)))
  vv <- c()
  for(i in 1:length(col_name))
  {
    vv <- c(vv, rep(col_name[i], 8))  #8 cell types
  }
  sample_col <- vv
  
  #3  
  celltype_col <- rep(c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts'), ncol(ictd_prop))
  #4
  #prediction_col <- rep(runif(8,min=0.01,max=0.8), ncol(ictd_prop))  
  prediction_col <- c()
  cell_type_order <- c('CD4.T.cells','CD8.T.cells','NK.cells','B.cells','monocytic.lineage','neutrophils','endothelial.cells','fibroblasts')
  for(i in 1:length(col_name))
  {
    prop_tmp <- ictd_prop[cell_type_order,i]
    prediction_col <- c(prediction_col, prop_tmp)
  }
  
  
  #bind
  output_df <- data.frame(dataset.name=dataset.name,sample.id=sample_col,cell.type=celltype_col,prediction=prediction_col)
  
  return(output_df)
}

multi_weight <- function(max_score)
{
  max_score_weight <- max_score
  #B cell
  max_score_weight[1] <- max_score_weight[1]*3
  max_score_weight[2] <- max_score_weight[2]*2
  max_score_weight[3] <- max_score_weight[3]*1.5
  #Fibro
  max_score_weight[5] <- max_score_weight[5]*3
  max_score_weight[6] <- max_score_weight[6]*2
  max_score_weight[7] <- max_score_weight[7]*1.5
  #Endo
  max_score_weight[9] <- max_score_weight[9]*3
  max_score_weight[10]<- max_score_weight[10]*2
  #Mono
  #max_score_weight[12]<- max_score_weight[12]*3
  #Neutro
  
  
  return(max_score_weight) 
}



pick_good_basis <- function(data.matrix, tg_R1_lists, loca_in_R4_tmp)
{
  max_cor_loca <- 1
  max_cor_value <- 0
  for(i in 1:length(loca_in_R4_tmp)){
    gn_tmp <- list()
    gn_tmp[[1]] <- tg_R1_lists[[loca_in_R4_tmp[i]]]
    names(gn_tmp) <- c("tmp_gene_list")
    svd_row_basis <- Compute_Rbase_SVD(data.matrix, gn_tmp)
    all_cor <- mean(cor(t(data.matrix[gn_tmp[[1]],]), t(svd_row_basis)))
    #print(paste('i=',i,' , average correlation = ',all_cor, sep = ""))
    if(all_cor == max_cor_value){ print("two candidate list have same svd-cor score!") }
    if(all_cor > max_cor_value){
      max_cor_value <- all_cor
      max_cor_loca <- i
    }
    
  }
  
  return(tg_R1_lists[[loca_in_R4_tmp[max_cor_loca]]])
}



select_R4_8_cellmk <- function(data.matrix,tg_R1_lists)
{
  load('TCGA_IM_core_markers.RData')
  #length(TCGA_core_list)
  names(TCGA_core_list)
  #length(TCGA_N_core_list)
  names(TCGA_N_core_list)
  print('Load TCGA marker done.')
  #1. B cell
  B_C_level_3 <- TCGA_core_list[[7]]
  B_C_level_2 <- c(TCGA_core_list[[18]], TCGA_core_list[[21]])
  B_N_level_1.5 <- TCGA_N_core_list[[14]]
  B_all <- c(B_C_level_3,B_C_level_2,B_N_level_1.5)
  #2. Fibroblast
  Fibro_C_level_3 <- TCGA_core_list[[4]]
  Fibro_C_level_2 <- c(TCGA_core_list[[9]],TCGA_core_list[[14]])
  Fibro_N_level_1.5 <- TCGA_N_core_list[[11]]
  Fibro_all <- c(Fibro_C_level_3,Fibro_C_level_2,Fibro_N_level_1.5)
  #3. Endothelial
  Endo_C_level_3 <- TCGA_core_list[[12]]
  Endo_N_level_2 <- TCGA_N_core_list[[21]]
  Endo_all <- c(Endo_C_level_3,Endo_N_level_2)
  #4. Monocyte
  Mono_C_level_3 <- TCGA_core_list[[5]]
  #5. Neutrophils
  Neutro_C_level_2 <- c(TCGA_core_list[[11]],TCGA_core_list[[13]],TCGA_core_list[[15]])
  Neutro_N_level_2 <- TCGA_N_core_list[[20]]
  #6.7.8 CD4T/CD8T/NK
  CD4T_C_level_3 <- TCGA_core_list[[8]]
  CD4T_N_level_2 <- TCGA_N_core_list[[3]]
  T_C_level_3 <- TCGA_core_list[[19]]
  T_N_level_3 <- TCGA_N_core_list[[19]]
  NK_CD8T_N_level_3 <- TCGA_N_core_list[[28]]
  
  TCGA_genegroup <- list(B_C_level_3,B_C_level_2,B_N_level_1.5,B_all,Fibro_C_level_3,Fibro_C_level_2,Fibro_N_level_1.5,Fibro_all,
                         Endo_C_level_3,Endo_N_level_2,Endo_all,Mono_C_level_3,Neutro_C_level_2,Neutro_N_level_2,
                         CD4T_C_level_3,CD4T_N_level_2,T_C_level_3,T_N_level_3,NK_CD8T_N_level_3)
  names(TCGA_genegroup) <- c("B_C_level_3","B_C_level_2","B_N_level_1.5","B_all","Fibro_C_level_3","Fibro_C_level_2","Fibro_N_level_1.5","Fibro_all",
                             "Endo_C_level_3","Endo_N_level_2","Endo_all","Mono_C_level_3","Neutro_C_level_2","Neutro_N_level_2",
                             "CD4T_C_level_3","CD4T_N_level_2","T_C_level_3","T_N_level_3","NK_CD8T_N_level_3")
  
  score_mat <- matrix(0, length(tg_R1_lists), 19)
  colnames(score_mat) <- names(TCGA_genegroup)
  rownames(score_mat) <- names(tg_R1_lists)
  for(i in 1:length(tg_R1_lists)){
    for(j in 1:length(TCGA_genegroup)){
      score_mat[i,j] <- length(intersect(tg_R1_lists[[i]],TCGA_genegroup[[j]])) / length(TCGA_genegroup[[j]])
    }
  }
  max_score <- matrix(apply(score_mat,2,max),1,ncol(score_mat))
  colnames(max_score) <- colnames(score_mat)
  #print(max_score)
  max_score_weight <- multi_weight(max_score)
  #print(max_score_weight)
  #use score to identify row ID which represent the true marker
  #1.B 
  B_loca <- which.max(max_score_weight[1:3])
  #B_mk <- tg_R1_lists[[which(score_mat[,B_loca] == max_score[,B_loca])[1]]]
  loca_in_R4_tmp <- which(score_mat[,B_loca] == max_score[,B_loca])
  B_mk <- pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp)
  #assume pick the first element
  # #----------how to pick if return two id?
  # tg_R1_lists[[85]]
  # [1] "CD79A" "MS4A1" "BANK1" "CD79B" "CD22"  "CD19"  "PAX5" 
  # > tg_R1_lists[[97]]
  # [1] "MS4A1"   "BANK1"   "CD79A"   "CD22"    "FAM129C" "FCRL1"   "FCRL5"  
  # #----------
  #2.Fibro
  Fibro_loca <- 5-1 + which.max(max_score_weight[5:8])
  loca_in_R4_tmp <- which(score_mat[,Fibro_loca] == max_score[,Fibro_loca])
  #Fibro_mk <- tg_R1_lists[[which(score_mat[,Fibro_loca] == max_score[,Fibro_loca])[1]]]
  Fibro_mk <- pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp)
  #assume pick the first element
  # > tg_R1_lists[[13]]
  # [1] "COL1A2" "COL1A1" "COL3A1" "LUM"    "THBS2"  "COL6A3" "MXRA8" 
  # > tg_R1_lists[[51]]
  # [1] "COL1A2" "COL1A1" "COL3A1" "PCOLCE" "THBS2"  "LUM"    "COL5A1"
  # > tg_R1_lists[[84]]
  # [1] "COL5A2" "COL1A2" "COL6A1" "COL6A3" "COL5A1" "ADAM12" "MMP2"  
  #3. Endo
  Endo_loca <- 9-1 + which.max(max_score_weight[9:11])
  #Endo_mk <- tg_R1_lists[[which(score_mat[,Endo_loca] == max_score[,Endo_loca])[1]]]
  loca_in_R4_tmp <- which(score_mat[,Endo_loca] == max_score[,Endo_loca])
  Endo_mk <- pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp)
  #4.Monocyte
  #Mono_mk <- tg_R1_lists[[which(score_mat[,12] == max_score[,12])[[1]]]]
  loca_in_R4_tmp <- which(score_mat[,12] == max_score[,12])
  Mono_mk <- pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp)
  # > tg_R1_lists[[19]]
  # [1] "MS4A6A" "CD14"   "AIF1"   "C1QA"   "C1QB"   "MS4A4A" "CD163" 
  # > tg_R1_lists[[68]]
  # [1] "CD14"    "CD163"   "MS4A6A"  "CSF1R"   "MSR1"    "MS4A4A"  "SIGLEC1"
  #5. Neutrophils
  Neutro_loca <- 13-1 + which.max(max_score_weight[13:14])
  #Neutro_mk <- tg_R1_lists[[which(score_mat[,Neutro_loca] == max_score[,Neutro_loca])[[1]]]]
  loca_in_R4_tmp <- which(score_mat[,Neutro_loca] == max_score[,Neutro_loca])
  Neutro_mk <- pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp)
  #6.7.8 T/NK
  TNK_mk <- c()
  for(k in 15:19){
    #mk_tmp <- tg_R1_lists[[which(score_mat[,k] == max_score[,k])[[1]]]]
    #TNK_mk <- c(TNK_mk, mk_tmp)
    loca_in_R4_tmp <- which(score_mat[,k] == max_score[,k])
    TNK_mk <- c(TNK_mk,pick_good_basis(data.matrix, tg_R1_lists, loca_in_R4_tmp))
  }
  
  eight_cell_mk <- list(B_mk,Fibro_mk,Endo_mk,Mono_mk,Neutro_mk,TNK_mk)
  names(eight_cell_mk) <- c("B.cells","fibroblasts","endothelial.cells","monocytic.lineage","neutrophils","TNK")
  
  print("select_R4_8_cellmk DONE!!!")
  return(eight_cell_mk)
}

# seperate_T4_8_NK <- function(data.matrix, eight_mk_list)
# {
#   #load LM22 signature
#   LM22 <- read.table("LM22.txt", sep="\t", header=T)
#   rownames(LM22) <- LM22[,1]
#   LM22 <- LM22[,-1]
#   LM22 <- as.matrix(LM22)
#   LM22_score <- cal_Zscore_signature(LM22)
#   
#   #1. CD4T
#   LM22_CD4T_col <- 5
#   LM22_CD4T_list1 <- names(which(LM22_score[,LM22_CD4T_col] == 1))
#   LM22_CD4T_col <- 6
#   LM22_CD4T_list2 <- names(which(LM22_score[,LM22_CD4T_col] == 1))
#   LM22_CD4T_col <- 7
#   LM22_CD4T_list3 <- names(which(LM22_score[,LM22_CD4T_col] == 1))
#   LM22_CD4T_all <- intersect(c(LM22_CD4T_list1,LM22_CD4T_list2,LM22_CD4T_list3), rownames(data.matrix))
#   LM22_CD4_tg_list <- list(LM22_CD4T_all)
#   names(LM22_CD4_tg_list) <- c("CD4T_all")
#   
#   tProp_CD4T <- Compute_Rbase_SVD(data.matrix, LM22_CD4_tg_list)
#   
#   cor(t(data.matrix[eight_cell_mk[[6]],]), t(tProp_CD4T) )
#   
#   eight_mk_list[[6]]
#   
#   monocyte_LM_loca <- 13
#   
# }

EPIC_NKT48_mk <- function()
{
  sigGeneEpic <- EPIC::TRef$sigGenes
  ref <- EPIC::TRef$refProfiles
  hhh <- ref[sigGeneEpic,]
  epic_cell_marker <- extract_marker(hhh,method="EPIC",add=T)
  
  TNK3_list <- list(epic_cell_marker[["CD4_mark"]], epic_cell_marker[["CD8_mark"]], epic_cell_marker[["NK_mark"]])
  names(TNK3_list) <- c("CD4.T.cells","CD8.T.cells","NK.cells")
  
  return(TNK3_list)
}

extract_combine_list <- function(mylist)
{
  comb <- c()
  for(i in 1:length(mylist))
  {
    comb <- c(comb, mylist[[i]])
  }
  
  return(comb)
}

filt_tg_list_by_dataName <- function(two_list_combine,data.matrix)
{
  dt_gn_name <- rownames(data.matrix)
  for(i in 1:length(two_list_combine))
  {
    two_list_combine[[i]] <- intersect(two_list_combine[[i]], dt_gn_name)
  }
  return(two_list_combine)
  
}

rm_NA_row <- function(data0)
{
  
  data1 <- data0[complete.cases(data0), ]
  print(paste("rm ", nrow(data0)-nrow(data1), " --rows", sep=""))
  
  return(data1)
}

normalize_data2_rm_sdZero <- function (data0) 
{
  data1 <- data0
  for (i in 1:nrow(data1)) {
    data1[i, ] <- (data1[i, ] - mean(data1[i, ]))/sd(data1[i, ])
  }
  data2 <- rm_NA_row(data1)
  return(data2)
}

BCV_ttest2_NA<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
{
  #print("call new BCV test2 function !!!")
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
  #print(ddd)
  
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


generate_unique_list<-function(tg_list,n_cut=1)
{
  aaa<-c()
  for(i in 1:length(tg_list))
  {
    aaa<-c(aaa,tg_list[[i]])
  }
  tg_selected<-names(which(table(aaa)<=n_cut))
  tg_list0<-tg_list
  for(i in 1:length(tg_list0))
  {
    tg_list0[[i]]<-intersect(tg_list[[i]],tg_selected)
  }
  return(tg_list0)
}

compute_jaccard2<-function(data_list1,data_list2)
{
  aaa<-matrix(0,length(data_list1),length(data_list2))
  rownames(aaa)<-names(data_list1)
  colnames(aaa)<-names(data_list2)
  for(i in 1:length(data_list1))
  {
    for(j in 1:length(data_list2))
    {
      
      aaa[i,j]<-length(intersect(data_list1[[i]],data_list2[[j]]))/min(length(data_list1[[i]]),length(data_list2[[j]]))
    }
  }
  return(aaa)
}

top5_SVD_cor<-function(data_c,tg_genes)
{
  tg_genes0<-intersect(rownames(data_c),tg_genes)
  tg_s2<-c()
  if(length(tg_genes0)>=5)
  {
    fff<-svd(data_c[tg_genes0,])$v[,1]
    if(mean(cor(fff,t(data_c[tg_genes0,])),ra.rm=T)<0)
    {
      fff<--fff
    }
    tg_s<-names(which(cor(fff,t(data_c[tg_genes0,]))[1,]>0.4))
    if(length(tg_s)>=5)
    {
      fff2<-svd(data_c[tg_s,])$v[,1]
      if(mean(cor(fff2,t(data_c[tg_s,])),ra.rm=T)<0)
      {
        fff2<--fff2
      }
      tg_s2<-names(sort(cor(fff2,t(data_c[tg_s,]))[1,],decreasing=T))[1:5]
    }
    else
    {
      tg_s2<-names(sort(cor(fff,t(data_c[tg_genes0,]))[1,],decreasing=T))[1:5]
    }
  }
  else
  {
    if(length(tg_genes0)>=3)
    {
      fff<-svd(data_c[tg_genes0,])$v[,1]
      if(mean(cor(fff,t(data_c[tg_genes0,])),ra.rm=T)<0)
      {
        fff<--fff
      }
      tg_s2<-names(which(cor(fff,t(data_c[tg_genes0,]))[1,]>0))
    }
  }
  return(tg_s2)
}





ICTD_round1 <- function(data_bulk)
{
  #1
  # load('TCGA-COAD_FPKM_T.RData')    #load COAD to test
  # rownames(data_t) <- gsub('^.*\\|', '', rownames(data_t))
  # rowname_stay <- setdiff(unique(rownames(data_t)), "" )
  # data_t <- data_t[rowname_stay,]
  # data_bulk <- data_t
  # #2
  # data_bulk <- GSE72056_diri_example[[1]]
  # 
  data.matrix = data_bulk
  if (length(colnames(data.matrix)) == 0) {
    warning("input data do NOT have colnames")
    colnames(data.matrix) <- paste("Setsample", 1:ncol(data.matrix),sep = "")
  }
  data.matrix <- rm_zero_row(data.matrix)
  if (max(data.matrix) > 30) {
    d.matrix <- log(data.matrix + 1)
    d.matrix <- as.matrix(d.matrix)
    print("do log because data > 20!")
  }else{
    d.matrix <- as.matrix(data.matrix)
  }
  data0 <- d.matrix
  data2 <- data0
  data21 <- data2[order(-apply(data2, 1, mean)), ]
  data22 <- data21[unique(rownames(data21)), ]
  data22_code <- data22[intersect(rownames(data22), TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3] == "protein_coding" & TCGA_ensem_annotation[,4] != ""), 4]), ]
  data23 <- data22[intersect(rownames(data22), TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3] == "protein_coding" & TCGA_ensem_annotation[,4] != ""), 4]), ]
  #data23 <- normalize_data2(data23)
  data23 <- normalize_data2_rm_sdZero(data23)
  data_CORS_cancer <- data23
  data_ccc <- data23
  list_c1 <- MRHCA_IM_compute_MR(data_CORS_cancer = data_ccc, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  MR_IM_result_new_c <- MRHCA_IM_compute_full_pub_new(data_CORS_cancer = data_ccc, list_c = list_c1, IM_id_list, immune_cell_uni_table = immune_cell_uni_table0_GS)
  tg_key = "nonono"
  
  list_new_c1<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key0,cell_type_enrich_cut=0.4,cor_cut0=0.8,num_cut=6,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
  R1_filter_step1_results_new1<-R1_list_filtering_step1_new_BCV2(list_new_c1,data_CORS_cancer=data_ccc,max_cut=10,cutn0=6,cut10=0.8,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)#cut10=0.8 for RNA-seq, cut10=0.75 for microarray, cut10=0.85 for scRNA-seq simulated data
  
  aaa1<-compute_jaccard2(R1_filter_step1_results_new1[[4]],LM22_test_list1_plus)
  
  list_new_c2<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key0,cell_type_enrich_cut=0.4,cor_cut0=0.75,num_cut=6,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
  R1_filter_step1_results_new2<-R1_list_filtering_step1_new_BCV2(list_new_c2,data_CORS_cancer=data_ccc,max_cut=10,cutn0=6,cut10=0.75,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)#cut10=0.8 for RNA-seq, cut10=0.75 for microarray, cut10=0.85 for scRNA-seq simulated data
  
  aaa2<-compute_jaccard2(R1_filter_step1_results_new2[[4]],LM22_test_list1_plus)
  
  list_new_c3<-Process_MR_IM_result_new(MR_IM_result_new_c,tg_key_c=tg_key0,cell_type_enrich_cut=0.25,cor_cut0=0.65,num_cut=6,num_cut2=5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)
  R1_filter_step1_results_new3<-R1_list_filtering_step1_new_BCV2(list_new_c3,data_CORS_cancer=data_ccc,max_cut=6,cutn0=6,cut10=0.6,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GS)#cut10=0.8 for RNA-seq, cut10=0.75 for microarray, cut10=0.85 for scRNA-seq simulated data
  
  aaa3<-compute_jaccard2(R1_filter_step1_results_new3[[4]],LM22_test_list1_plus)
  
  #M1----------------------------------------------------------------------------------------
  tg_ids<-c(1:8)
  tg_markers<-list()
  for(i in 1:length(tg_ids))
  {
    #tg_CCC<-list()
    st<-0
    if(sum(aaa1[,tg_ids[i]]>0.7)>0)
    {
      tg_id_c<-which(aaa1[,tg_ids[i]]>0.7)
      if(length(tg_id_c)==1)
      {
        tg_markers[[i]]<-R1_filter_step1_results_new1[[4]][[tg_id_c]]
        st<-1
      }
      if(length(tg_id_c)>1)
      {
        nn<-c()
        for(j in 1:length(tg_id_c))
        {
          nn<-c(nn,length(intersect(LM22_test_list1_plus_unique_core[[tg_ids[i]]],R1_filter_step1_results_new1[[4]][[tg_id_c[j]]])))
        }
        if(max(nn)>0)
        {
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]####################################
          tg_markers[[i]]<-R1_filter_step1_results_new1[[4]][[tg_id_ccc]]
          st<-1
        }
        else
        {
          nn<-c()
          for(j in 1:length(tg_id_c))
          {
            nn<-c(nn,mean(cor(t(data_CORS_cancer[R1_filter_step1_results_new1[[4]][[tg_id_c[j]]],]))))
          }
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]#####################################
          tg_markers[[i]]<-R1_filter_step1_results_new1[[4]][[tg_id_ccc]]
          st<-1
        }
      }
    }
    if((st==0)&(sum(aaa2[,tg_ids[i]]>0.7)>0))
    {
      tg_id_c<-which(aaa2[,tg_ids[i]]>0.7)
      if(length(tg_id_c)==1)
      {
        tg_markers[[i]]<-R1_filter_step1_results_new2[[4]][[tg_id_c]]
        st<-1
      }
      if(length(tg_id_c)>1)
      {
        nn<-c()
        for(j in 1:length(tg_id_c))
        {
          nn<-c(nn,length(intersect(LM22_test_list1_plus_unique_core[[tg_ids[i]]],R1_filter_step1_results_new2[[4]][[tg_id_c[j]]])))
        }
        if(max(nn)>0)
        {
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]##################################
          tg_markers[[i]]<-R1_filter_step1_results_new2[[4]][[tg_id_ccc]]
          st<-1
        }
        else
        {
          nn<-c()
          for(j in 1:length(tg_id_c))
          {
            nn<-c(nn,mean(cor(t(data_CORS_cancer[R1_filter_step1_results_new2[[4]][[tg_id_c[j]]],]))))
          }
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]#####################################
          tg_markers[[i]]<-R1_filter_step1_results_new2[[4]][[tg_id_ccc]]
          st<-1
        }
      }
    }
    if((st==0)&(sum(aaa3[,tg_ids[i]]>0.7)>0))
    {
      tg_id_c<-which(aaa3[,tg_ids[i]]>0.7)
      if(length(tg_id_c)==1)
      {
        tg_markers[[i]]<-R1_filter_step1_results_new3[[4]][[tg_id_c]]
        st<-1
      }
      if(length(tg_id_c)>1)
      {
        nn<-c()
        for(j in 1:length(tg_id_c))
        {
          nn<-c(nn,length(intersect(LM22_test_list1_plus_unique_core[[tg_ids[i]]],R1_filter_step1_results_new3[[4]][[tg_id_c[j]]])))
        }
        if(max(nn)>0)
        {
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]###########################################
          tg_markers[[i]]<-R1_filter_step1_results_new3[[4]][[tg_id_ccc]]
          st<-1
        }
        else
        {
          nn<-c()
          for(j in 1:length(tg_id_c))
          {
            nn<-c(nn,mean(cor(t(data_CORS_cancer[R1_filter_step1_results_new3[[4]][[tg_id_c[j]]],]))))
          }
          tg_id_ccc<-tg_id_c[which(nn==max(nn))[1]]#########################################
          tg_markers[[i]]<-R1_filter_step1_results_new3[[4]][[tg_id_ccc]]
          st<-1
        }
      }
    }
    if(st==0)
    {
      tg_markers[[i]]<-0
    }
  }
  names(tg_markers)<-names(LM22_test_list1_plus_unique_core)[tg_ids]
  for(i in 1:length(tg_markers))
  {
    if(length(tg_markers[[i]])==1)
    {
      tg_markers[[i]]<-top5_SVD_cor(data01,LM22_test_list1_plus_unique_core[[tg_ids[i]]]) ##################################################
    }
  }
  
  
  l4<-tg_markers
  b4<-Compute_Rbase_SVD(data.matrix,l4)
  colnames(b4)<-colnames(data.matrix)
  Prop <- b4
  rownames(Prop) <- c('B.cells','CD4.T.cells','CD8.T.cells','NK.cells','neutrophils','monocytic.lineage','fibroblasts','endothelial.cells')
  
  print("ICTD_round1 DONE!!!")
  return(Prop)
  
  
  #NMF  
  #NMF test and debug
  #how to build matrix C?
  #test and debug
  
  
}

#----------------------function part finish---------




#--------------pipeline-------------
print("test ictd_round1!")


print(list.files())
print('current dir:')
print(getwd())
print("what in the input folder:")
print(list.files('input/'))
input_df <- read.csv('input/input.csv')
print(input_df)
print('read <input.csv> done!')

load("LM22_marker_list.RData")
print('load LM22 list done!')

# Extract the names of each dataset
dataset_name <- as.character(input_df$dataset.name)

# Extract the names of the expression files that use gene name
expression_files <- as.character(input_df$hugo.expr.file)

input_combine <- c()
output_all_ds <- c()
for(i in 1:length(expression_files))
{
  ff_tmp <- paste('input/', expression_files[i],sep='')
  print(ff_tmp)
  data_tmp <- read.csv(ff_tmp)
  rownames(data_tmp) <- data_tmp[,1]
  data_tmp <- data_tmp[,-1]
  #data_bulk <- data_tmp
  data_tmp <- as.matrix(data_tmp)
  data_tmp[is.na(data_tmp)] <- 0
  
  print(paste(ff_tmp, " dim is :"))
  print(dim(data_tmp))
  print(data_tmp[1:5,1:5])
  ictd_result <- ICTD_round1(data_tmp)
  ictd_prop <- ictd_result
  
  
  # # a)fake output!!!
  # ictd_prop <- data_tmp[1:12,]
  # dn_tmp <- unlist(strsplit(expression_files[i],split='.',fixed=T))[[1]]
  # output_tmp <- ictd_2_output(ictd_prop, dn_tmp)
  # #combine prediction into big dataframe
  # output_all_ds <- rbind(output_all_ds, output_tmp)
  #b) real output!!!
  dn_tmp <- unlist(strsplit(expression_files[i],split='.',fixed=T))[[1]]
  output_tmp <- ictd_2_output_real(ictd_prop, dn_tmp)
  #combine prediction into big dataframe
  output_all_ds <- rbind(output_all_ds, output_tmp)
}

# Create the directory the output will go into
dir.create("output")

# Write the result into output directory
readr::write_csv(output_all_ds, "output/predictions.csv")
print("files exist in the output folder:")
print(list.files("output"))
print("output file dim:")
print(dim(output_all_ds))
print("output :::")
print(output_all_ds)



