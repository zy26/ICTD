#Functions
###################
#filter_1rank_genes<-function(tg_cluster_genes,data_CORS_cancer,cut0=0.8)
#compute_adjusted_cor_all<-function(tg_cluster_genes,data_CORS_cancer,qanth=0.7,qantl=0.3)
#compute_adjusted_cor_higher<-function(tg_cluster_genes,data_CORS_cancer,qant=0.7)
#compute_adjusted_cor_lower<-function(tg_cluster_genes,data_CORS_cancer,qant=0.3)
#compute_BCV_1fold_test<-function(marker_list,data_c,ccc0)
#compute_cor<-function(data_list)
#compute_cor2<-function(data_list,data_t)
#compute_base2<-function(data_list,data_t)
#compute_min_jaccard<-function(data_list)
#BCV_ttest<-function(data0,rounds=20,slice0=2,maxrank0=30)
#evaluate_rank1_markers<-function(data_c,tg_marker_lists)
#compute_O_Rowspace<-function(data0,k)
#compute_O_Rowspace0<-function(data0,k)
#MRHCA_IM_compute_MR<-function(data_CORS_cancer,IM_id_list,immune_cell_uni_table=marker_stats20_uni,step_size0=50)
#MRHCA_IM_compute_full_pub<-function(data_CORS_cancer,list_c,IM_id_list,immune_cell_uni_table=marker_stats20_uni,step_size0=20)
#plot_simple1<-function(tg_gene,MR_M,cor_c,tg_id,ll=200,immune_cell_uni_table=immune_cell_uni_table0)
#scoring_MR_order<-function(aaa,penalty=-2,lbb=-20)
#calculate_growth_rate2<-function(x0,step=20)
#Process_MR_IM_result<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5,num_cut=5)
#Process_MR_IM_result_GPL570<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5)
#BCV_selection_top_bases<-function(data_CORS_cancer,tg_key_c=tg_key_c,tg_1_rank_markers=tg_1_rank_markers,immune_cell_uni_table=immune_cell_uni_table0,hcutn=40)
#BCV_screen_top_bases<-function(data_CORS_cancer,tg_key_c=tg_key_c,tg_1_rank_markers=tg_1_rank_markers,immune_cell_uni_table=immune_cell_uni_table0,hcutn=40)
#hcut_update<-function(tg_trunc_list,tg_trunc_list_old)
#plot_simple2<-function(MR_list,tg_gene="",ll=200)
#MAP_GPL570_genes2<-function(tg_cluster_genes)
#cut_modules<-function(tg_list,cutn=20)
#MAP_GPL570_genes<-function(tg_cluster_genes)
#R1_list_filtering_step1<-function(list_c2,data_CORS_cancer,data01,data0,cutn0=15,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
#R1_list_filtering_step2_BCV<-function(R1_marker_list_f2,data_CORS_cancer)
#clean_rank1_module<-function(data_c,module_info,module_rank,st0=8)
#R1_list_filtering_step3_buildNMF<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
#BCV_1fold_test_fast<-function(marker_list,data_c)
#compute_CompRowspace_NN<-function(data_cc=tg_data,module_RB=module_selected[[4]],ROUNDS=3)
#Select_cancer_module<-function(cancer_module,data_cc=data_t,kk=1)
#BCV_1fold_test_fast_diag<-function(marker_list,data_c)
#R1_list_filtering_step25_build_LG<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
#R1_list_filtering_step35_buildNMF_pure1<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
#R1_list_filtering_step35_buildNMF_pure2<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
#generate_oc_1c_cor<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
#R1_list_filtering_step3_buildNMF_new<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
#Reduce_NMF_redundancy<-function(C_ids,NMF_indi_all,data_CORS_cancer)
#BCV_ttest2<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
#merge_R1_2<-function(RBases_R1_2,R1_list,cor_cut=0.9)

###############
filter_1rank_genes<-function(tg_cluster_genes,data_CORS_cancer,cut0=0.8)
{
	tg_list_new<-list()
	for(i in 1:length(tg_cluster_genes))
	{
		nn<-min(length(tg_cluster_genes[[i]]),10)
		ccc<-RMSE_row(data_CORS_cancer[tg_cluster_genes[[i]],]%*%compute_O_Rowspace0(data_CORS_cancer[tg_cluster_genes[[i]][1:nn],],1))/RMSE_row(data_CORS_cancer[tg_cluster_genes[[i]],])
		tg_list_new[[i]]<-names(which(ccc>cut0))
		if(length(tg_list_new[[i]])<3)
		{
			tg_list_new[[i]]<-unique(c(tg_list_new[[i]],tg_cluster_genes[[i]][1:3]))[1:3]
		}
	}
	names(tg_list_new)<-names(tg_cluster_genes)
	return(tg_list_new)
}
	compute_adjusted_cor_all<-function(tg_cluster_genes,data_CORS_cancer,qanth=0.7,qantl=0.3)
	{
		tg_cluster_bases<-compute_base2(tg_cluster_genes,data_CORS_cancer)
		cor_all<-c()
		for(i in 1:length(tg_cluster_genes))
		{
			tg_samples<-which(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean)>quantile(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean),qanth))
			cor_c<-cor(t(tg_cluster_bases[,tg_samples]))
			cor_all<-rbind(cor_all,as.vector(cor_c))
		}
		for(i in 1:length(tg_cluster_genes))
		{
			tg_samples<-which(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean)<quantile(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean),qantl))
			cor_c<-cor(t(tg_cluster_bases[,tg_samples]))
			cor_all<-rbind(cor_all,as.vector(cor_c))
		}
		cor_0<-apply(cor_all,2,min)
		cor_1<-matrix(as.vector(cor_0),length(tg_cluster_genes),length(tg_cluster_genes))
		return(cor_1)
	}

	compute_adjusted_cor_higher<-function(tg_cluster_genes,data_CORS_cancer,qant=0.7)
	{
		tg_cluster_bases<-compute_base2(tg_cluster_genes,data_CORS_cancer)
		cor_all<-c()
		for(i in 1:length(tg_cluster_genes))
		{
			tg_samples<-which(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean)>quantile(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean),qant))
			cor_c<-cor(t(tg_cluster_bases[,tg_samples]))
			cor_all<-rbind(cor_all,as.vector(cor_c))
		}
		cor_0<-apply(cor_all,2,min)
		cor_1<-matrix(as.vector(cor_0),length(tg_cluster_genes),length(tg_cluster_genes))
		return(cor_1)
	}
	
	compute_adjusted_cor_lower<-function(tg_cluster_genes,data_CORS_cancer,qant=0.3)
	{
		tg_cluster_bases<-compute_base2(tg_cluster_genes,data_CORS_cancer)
		cor_all<-c()
		for(i in 1:length(tg_cluster_genes))
		{
			tg_samples<-which(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean)<quantile(apply(data_CORS_cancer[tg_cluster_genes[[i]],],2,mean),qant))
			cor_c<-cor(t(tg_cluster_bases[,tg_samples]))
			cor_all<-rbind(cor_all,as.vector(cor_c))
		}
		cor_0<-apply(cor_all,2,min)
		cor_1<-matrix(as.vector(cor_0),length(tg_cluster_genes),length(tg_cluster_genes))
		return(cor_1)
	}



compute_BCV_1fold_test<-function(marker_list,data_c,ccc0)
{
	stat<-matrix(0,length(marker_list),length(marker_list))
	rownames(stat)<-names(marker_list)
	colnames(stat)<-names(marker_list)
	stat_i<-stat
	stat_j<-stat
	stat_k<-stat
	i<-"ITGAL"
	j<-"C1QB"
	print("Merge 1Rank Markers BCV test!")
	for(i in 1:length(marker_list))
	{
		print(c(i,length(marker_list)))
		for(j in 1:length(marker_list))
		{
			if(i<j)
			{	
				ccc1<-data_c[marker_list[[i]],]
				ccc2<-data_c[marker_list[[j]],]
				ccc<-rbind(ccc1,ccc2)
				ccc<-ccc[unique(rownames(ccc)),]
				tg_samples1<-which(apply(ccc1,2,mean)>quantile(apply(ccc1,2,mean),0.5))
				tg_samples2<-which(apply(ccc2,2,mean)>quantile(apply(ccc2,2,mean),0.5))
				i1<-sum(BCV_ttest(ccc1,maxrank0=10)<0.001)
				i2<-sum(BCV_ttest(ccc2,maxrank0=10)<0.001)
				i12<-sum(BCV_ttest(ccc,maxrank0=10)<0.001)
				j1<-sum(BCV_ttest(ccc1[,tg_samples1],maxrank0=10)<0.001)
				j2<-sum(BCV_ttest(ccc2[,tg_samples1],maxrank0=10)<0.001)
				j12<-sum(BCV_ttest(ccc[,tg_samples1],maxrank0=10)<0.001)
				k1<-sum(BCV_ttest(ccc1[,tg_samples2],maxrank0=10)<0.001)
				k2<-sum(BCV_ttest(ccc2[,tg_samples2],maxrank0=10)<0.001)
				k12<-sum(BCV_ttest(ccc[,tg_samples2],maxrank0=10)<0.001)
				stat_i[i,j]<-stat_i[j,i]<-i12-i1-i2
				stat_j[i,j]<-stat_j[j,i]<-j12-j1-j2
				stat_k[i,j]<-stat_k[j,i]<-k12-k1-k2
			}
		}
	}
	print("Merge 1Rank Markers BCV test done!")
	tg_link_table<-sign((stat_k<0)+(stat_j<0)+ccc0)
	tg_link_table[which(is.na(tg_link_table))]<-0
	hh<-heatmap.2(tg_link_table,Rowv=T,Colv =T,scale="none",main="",
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)
	tgs_or<-rownames(tg_link_table)[hh$rowInd]
tgs_o<-tgs_or[length(tgs_or):1]
order_list<-c()
for(i in 1:length(tgs_o))
{
	tg_ids_c<-tgs_o[i]
	for(j in i:length(tgs_o))
	{
		tg_ids_c10<-unique(c(tg_ids_c,tgs_o[j]))
		if(sum(tg_link_table[tg_ids_c10,tg_ids_c10])==(length(tg_ids_c10))^2)
		{
			tg_ids_c<-tg_ids_c10
		}
	}
	order_list[[i]]<-tg_ids_c
}
for(i in length(order_list):2)
{
	for(j in (i-1):1)
	{
		if(length(intersect(order_list[[i]],order_list[[j]]))==length(order_list[[i]]))
		{
			order_list[[i]]<-""
		}
	}
}
ccc<-c()
for(i in 1:length(order_list))
{
	if(sum(order_list[[i]]=="")==0)
	{
		ccc<-c(ccc,order_list[[i]])
	}
}
tg_genes<-names(which(table(ccc)>=4))
Rank_1_merged_modules<-list()
Rank_1_merged_hubs<-list()
N<-0
if(length(tg_genes)>0)
{
for(i in 1:length(tg_genes))
{
	N<-N+1
	Rank_1_merged_modules[[N]]<-marker_list[[tg_genes[i]]]
	Rank_1_merged_hubs[[N]]<-tg_genes[i]
}
}
for(i in 1:length(order_list))
{
	tg_ccc<-setdiff(order_list[[i]],c(tg_genes,""))
	if(length(tg_ccc)>0)
	{
		tg_module<-c()
		for(j in 1:length(tg_ccc))
		{
			tg_module<-unique(c(tg_module,marker_list[[tg_ccc[j]]]))
		}
		st<-0
		if(length(Rank_1_merged_modules)>0)
	      {
		for(j in 1:length(Rank_1_merged_modules))
		{
			if(length(intersect(Rank_1_merged_modules[[j]],tg_module))==length(tg_module))
			{
				st<-1
			}
		}
            }
		if(st==0)
		{
			N<-N+1
			Rank_1_merged_modules[[N]]<-tg_module
			Rank_1_merged_hubs[[N]]<-tg_ccc
		}
	}
}
	return(list(Rank_1_merged_hubs,Rank_1_merged_modules,stat_i,stat_k,stat_j))
}

compute_cor<-function(data_list)
{
	ccc<-c()
	for(i in 1:length(data_list))
	{
		ccc<-rbind(ccc,data_list[[i]])
	}
	return(cor(t(ccc)))
}

compute_cor2<-function(data_list,data_t)
{
	ccc<-c()
	for(i in 1:length(data_list))
	{
		ccc<-rbind(ccc,apply(normalize_data2(data_t[data_list[[i]],]),2,mean))
	}
	return(cor(t(ccc)))
}



compute_cor2<-function(data_list,data_t)
{
	ccc<-c()
	for(i in 1:length(data_list))
	{
		ccc<-rbind(ccc,apply(normalize_data2(data_t[data_list[[i]],]),2,mean))
	}
	return(cor(t(ccc)))
}



compute_base2<-function(data_list,data_t)
{
        ccc<-c()
        for(i in 1:length(data_list))
        {
                SVD_ccc<-svd(data_t[data_list[[i]],])
		    ddd<- t(SVD_ccc$v)[1,]

		    if(cor(ddd,apply(data_t[data_list[[i]],],2,mean))<0)
			{
				ddd<--ddd
			}	
               ccc<-rbind(ccc,ddd)
        }
        return(ccc)
}

compute_min_jaccard<-function(data_list)
{
        aaa<-matrix(0,length(data_list),length(data_list))
        rownames(aaa)<-names(data_list)
        colnames(aaa)<-names(data_list)
        for(i in 1:length(data_list))
        {
                for(j in 1:length(data_list))
                {
                        if(i<j)
                        {
                                aaa[i,j]<-aaa[j,i]<-length(intersect(data_list[[i]],data_list[[j]]))/min(length(data_list[[i]]),length(data_list[[j]]))
                        }
                }
        }
        diag(aaa)<-1
        return(aaa)
}

###########
BCV_ttest<-function(data0,rounds=20,slice0=2,maxrank0=30)
{
	  x<-data0
        fff_cc<-c()
        for(kk in 1:rounds)
        {
                cv_result <- cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
                fff_cc<-rbind(fff_cc,cv_result$msep)
        }
        pp<-c()
        for(kk in 1:(ncol(fff_cc)-1))
        {
                pp_c<-1
                if(mean(fff_cc[,kk],na.rm=T)>mean(fff_cc[,kk+1],na.rm=T))
                {
                        pp_c<-t.test(fff_cc[,kk],fff_cc[,kk+1])$p.value
                }
                pp<-c(pp,pp_c)
        }
	return(pp)       
}

evaluate_rank1_markers<-function(data_c,tg_marker_lists)
{
	tg_data<-data_c
	ddd<-c()
	rowsp1<-list()
	for(i in 1:length(tg_marker_lists))
	{
		tg_genes<-tg_marker_lists[[i]]
		nn<-length(tg_genes)
		ave_p<-(sum(cor(t(tg_data[tg_genes,])))-nn)/(((nn-1)*(nn)))
		data_cc<-tg_data[tg_genes,]
		fff<-t(svd(data_cc)$v)[1,]
		if(mean(cor(fff,t(data_cc)))<0)
		{
			fff<--fff
		}
		names(fff)<-colnames(data_cc)
		rowsp1[[i]]<-fff
		pp_c<-BCV_ttest(tg_data[tg_genes,])
		dim_stat_c<-length(which(pp_c<(0.001)))
		data_C_cc<-data_cc%*%compute_O_Rowspace0(data_cc,1)
		ccc1<-RMSE_row(data_C_cc)
		ccc2<-RMSE_row(data_cc)
		data_cc0<-normalize_data2(data_cc)
		data_C_cc0<-data_cc0%*%compute_O_Rowspace0(data_cc0,1)
		ccc10<-RMSE_row(data_C_cc0)
		ccc20<-RMSE_row(data_cc0)
		RMSE_t<-mean(ccc1/ccc2)
		RMSE_t_21<-mean(ccc1)/mean(ccc2)
		RMSE_adj<-mean(ccc10/ccc20)
		ccc<-c(ave_p,dim_stat_c,RMSE_t,RMSE_t_21,RMSE_adj)
		ddd<-rbind(ddd,ccc)
	}
	rownames(ddd)<-names(tg_marker_lists)
	names(rowsp1)<-names(tg_marker_lists)
	colnames(ddd)<-c("Ave_p","Dim","RMSE_t","RMSE_t_21","RMSE_adj")
	eee<-list(ddd,rowsp1)
	names(eee)<-c("Stat","RowBase")
	return(eee)
}	


compute_O_Rowspace<-function(data0,k)
{
        SVD_ccc<-svd(normalize_data2(data0))
        ttt_ccc<-t(SVD_ccc$v)[1,]%*%t(t(SVD_ccc$v)[1,])/sum((t(SVD_ccc$v)[1,])^2)
        t_ccc0<-ttt_ccc*0
        for(i in 1:k)
        {
                t_ccc0<-t_ccc0+t(SVD_ccc$v)[i,]%*%t(t(SVD_ccc$v)[i,])/sum((t(SVD_ccc$v)[i,])^2)
        }
        return(t_ccc0)
}

compute_O_Rowspace0<-function(data0,k)
{
        SVD_ccc<-svd(data0)
        ttt_ccc<-t(SVD_ccc$v)[1,]%*%t(t(SVD_ccc$v)[1,])/sum((t(SVD_ccc$v)[1,])^2)
        t_ccc0<-ttt_ccc*0
        for(i in 1:k)
        {
                t_ccc0<-t_ccc0+t(SVD_ccc$v)[i,]%*%t(t(SVD_ccc$v)[i,])/sum((t(SVD_ccc$v)[i,])^2)
        }
        return(t_ccc0)
}
####################

###########
MRHCA_IM_compute_MR<-function(data_CORS_cancer,IM_id_list,immune_cell_uni_table=marker_stats20_uni,step_size0=50)
{
	immune_cell_markers<-unique(rownames(immune_cell_uni_table))
	tg_genes_all<-intersect(rownames(data_CORS_cancer),immune_cell_markers)
	cor_c<-cor(t(data_CORS_cancer[tg_genes_all,]))
      diag(cor_c)<-0
      cor_r_rank1<-apply(cor_c,1,rank)
      cor_r_rank2<-t(apply(cor_c,2,rank))
      cor_r_rank1<-nrow(cor_r_rank1)-cor_r_rank1+1
      cor_r_rank2<-nrow(cor_r_rank2)-cor_r_rank2+1
      MR_M<-cor_r_rank1*cor_r_rank2
    	ddd<-Empirical_null_distribution_MR(data_CORS_cancer[tg_genes_all,],Rounds=1000)
      null_growth_rate<-c()
      step_size<-step_size0   #parameter
      print("Generating Empirical Null Growth Rate")
      for(i in 1:1000)
      {
             x0<-sqrt(sort(ddd[i,]))
             null_growth_rate<-rbind(null_growth_rate,calculate_growth_rate2(x0,step=step_size))
      }
      up_lim<-mean(apply(null_growth_rate,2,quantile,probs =0.995)[1:200])
      down_lim<-mean(apply(null_growth_rate,2,quantile,probs =0.005)[1:200])
	return(list(MR_M,cor_c,down_lim))
}

MRHCA_IM_compute_full_pub<-function(data_CORS_cancer,list_c,IM_id_list,immune_cell_uni_table=marker_stats20_uni,step_size0=20)
{
	MR_M<-list_c[[1]]
	cor_c<-list_c[[2]]
	down_lim<-list_c[[3]]
 	print("Compute MR IM genes")
     	MR_IM_result_c<-list()
	print(nrow(MR_M))
	step_size<-step_size0
	IM_id_list0<-IM_id_list
	MR_IM_result_c<-c()
      for(ii in 1:nrow(list_c[[1]]))
      {
		tg_gene<-rownames(MR_M)[ii]
		x0<-sqrt(sort(MR_M[tg_gene,]))
		tg_growth_rate<-calculate_growth_rate2(x0,step=step_size)
		tg_ccc1<-which(tg_growth_rate<down_lim)
		aaa<-sqrt(sort(MR_M[tg_gene,]))
		bbb3<-aaa/(1:length(aaa))
		bbb4<-cor_c[tg_gene,names(sort(MR_M[tg_gene,]))]
		ccc<-immune_cell_uni_table[names(sort(MR_M[tg_gene,])),]^(1/2)
		ccc0<-c()
		for(j in 1:length(IM_id_list0))
		{	
			if(length(IM_id_list[[j]])>1)
			{
				xx<-sum((1/(1:length(IM_id_list[[j]])))^(1/2))
				cc0<-(cumsum(apply(ccc[,IM_id_list[[j]]],1,sum))/xx/c(1:nrow(ccc)))
			}
			else
			{
				cc0<-(cumsum(ccc[,IM_id_list[[j]]])/c(1:nrow(ccc)))
			}
			ccc0<-cbind(ccc0,cc0)
		}
		colnames(ccc0)<-names(IM_id_list)
		ddd<-cbind(1:length(tg_growth_rate),tg_growth_rate,bbb3,bbb4,ccc0)
		colnames(ddd)[1:4]<-c("MR_order_ID","growth_rate","sorted_MR","Corr") 
		tg_ccc3<-intersect(tg_ccc1,which(bbb4>0.6))
  		fff<-""
		if(length(tg_ccc3)>1)
		{
			fff<-as.matrix(ddd[1:max(tg_ccc3),])
		}
		if(ii%%500==1)
		{
			print(ii)
		}
		MR_IM_result_c[[ii]]<-list(fff,tg_ccc3)
      }
      names(MR_IM_result_c)<-rownames(MR_M)[1:length(MR_IM_result_c)]
      return(MR_IM_result_c)
}

plot_simple1<-function(tg_gene,MR_M,cor_c,tg_id,ll=200,immune_cell_uni_table=immune_cell_uni_table0)
{
          aaa<-sqrt(sort(MR_M[tg_gene,]))#-(1:nrow(MR_M)) 
          bbb<-aaa/(1:length(aaa))
          plot(bbb,type="b",pch=16,xlim=c(1,ll),ylim=c(0,2),lwd=2,col="lightblue",main=tg_gene)
          bbb1<-bbb

          ccc<-immune_cell_uni_table[names(sort(MR_M[tg_gene,])),]^(1/2)
          points((cumsum(ccc[,tg_id])/c(1:nrow(ccc))),pch=16,type="b",lwd=2,col="palegreen")
          bbb2<-(cumsum(ccc[,tg_id])/c(1:nrow(ccc)))
 
          x0<-sqrt(sort(MR_M[tg_gene,]))
          tg_growth_rate<-calculate_growth_rate2(x0,step=200)
          points(tg_growth_rate,pch=16,type="b",lwd=2,col="violet",ylim=c(-1,2))
          bbb3<-tg_growth_rate

          ccc<-cor_c[tg_gene,names(sort(MR_M[tg_gene,]))]
          points(ccc,pch=16,type="b",lwd=2,col="orange")
          bbb4<-ccc
}



Process_MR_IM_result<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5,num_cut=5)
{
tg_1_rank_markers<-list()
tg_m_names<-c()
N<-0
print("Select Marker!")
for(i in 1:length(MR_IM_result_c))
{
        if(length(MR_IM_result_c[[i]][[2]])>=num_cut)
        {
                ss<-scoring_MR_order(MR_IM_result_c[[i]][[2]])
                if(ss>=num_cut)
                {
                        ccc<-MR_IM_result_c[[i]][[1]][1:ss,]
                        if(sum(ccc[,4]>cor_cut0)>=num_cut)
                        {       
                                tg_ccc<-names(which(ccc[,4]>cor_cut0))
                                if(length(tg_ccc)>=num_cut)
                                {
                                        ddd<-apply(ccc[,-c(1:4)],1,max)
                                        fff<-which(ddd>cell_type_enrich_cut)
                                        if(length(intersect(which(ccc[,4]>cor_cut0),fff))>=num_cut)
                                        {
                                                N<-N+1
                                                tg_1_rank_markers[[N]]<-ccc[which(ccc[,4]>cor_cut0),]
                                                tg_m_names<-c(tg_m_names,names(MR_IM_result_c)[i])
                                        }
                                }
                        }
                }
        }
}
print("Select Marker Done!")
#tg_RF2<-paste(tg_key_c,"_1rankmarker_cell_type_consistency.pdf",sep="")
#library(gplots)
#colors = c(0:100)/100
#my_palette <- grDevices::colorRampPalette(c("white","white", "blue"))(n =100)
#pdf(tg_RF2)    
#for(i in 1:length(tg_1_rank_markers))
#{
#       aaa<-tg_1_rank_markers[[i]][,-c(1:4)]
#       heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
#       col=my_palette,breaks=colors,density.info="none",dendrogram="both",
#       trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)
#}
#dev.off()
names(tg_1_rank_markers)<-tg_m_names
list_cc<-list(tg_1_rank_markers,tg_m_names)
return(list_cc)
}

Process_MR_IM_result_GPL570<-function(MR_IM_result_c=MR_IM_result_c,tg_key_c=tg_key_c,cor_cut0=0.7,cell_type_enrich_cut=0.5,num_cut=5)
{
tg_1_rank_markers<-list()
tg_m_names<-c()
N<-0
print("Select Marker!")
for(i in 1:length(MR_IM_result_c))
{
	if(length(MR_IM_result_c[[i]][[2]])>=num_cut)
	{
		ss<-scoring_MR_order(MR_IM_result_c[[i]][[2]])
		if(ss>=num_cut)
		{
			ccc<-MR_IM_result_c[[i]][[1]][1:ss,]
			if(sum(ccc[,4]>cor_cut0)>=num_cut)
			{	
				tg_ccc<-unique(GPL570_id_symbol0[intersect(GPL570_id_symbol[,1],c(names(MR_IM_result_c)[i],rownames(ccc[which(ccc[,4]>cor_cut0),]))),2])
				if(length(tg_ccc)>=num_cut)
				{
					ddd<-apply(ccc[,-c(1:4)],1,max)
					fff<-which(ddd>cell_type_enrich_cut)
					if(length(intersect(which(ccc[,4]>cor_cut0),fff))>=num_cut)
					{
						N<-N+1
						tg_1_rank_markers[[N]]<-ccc[which(ccc[,4]>cor_cut0),]
						tg_m_names<-c(tg_m_names,names(MR_IM_result_c)[i])
					}
				}
			}
		}
	}
}
print("Select Marker Done!")
#tg_RF2<-paste(tg_key_c,"_1rankmarker_cell_type_consistency.pdf",sep="")
#library(gplots)
#colors = c(0:100)/100
#my_palette <- grDevices::colorRampPalette(c("white","white", "blue"))(n =100)
#pdf(tg_RF2)	
#for(i in 1:length(tg_1_rank_markers))
#{
#	aaa<-tg_1_rank_markers[[i]][,-c(1:4)]
#	heatmap.2(aaa,Rowv=F,Colv =F,scale="none",main=tg_m_names[i],
#	col=my_palette,breaks=colors,density.info="none",dendrogram="both",
#	trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)
#}
#dev.off()
names(tg_1_rank_markers)<-tg_m_names
list_cc<-list(tg_1_rank_markers,tg_m_names)
return(list_cc)
}

BCV_selection_top_bases<-function(data_CORS_cancer,tg_key_c=tg_key_c,tg_1_rank_markers=tg_1_rank_markers,immune_cell_uni_table=immune_cell_uni_table0,hcutn=40)
{
library(bcv)
print("start BCV")
pp_all<-list()
cc_all<-c()
for(j in 1:length(tg_1_rank_markers))
{
	aaa<-data_CORS_cancer[c(names(tg_1_rank_markers)[j],rownames(tg_1_rank_markers[[j]])),]
	x<-aaa
	fff_cc<-c()
	for(kk in 1:10)
	{
		cv_result <- cv.svd.gabriel(x, 3, 3, maxrank = 30)
		fff_cc<-rbind(fff_cc,cv_result$msep)
	}
	pp<-c()
	for(kk in 1:(ncol(fff_cc)-1))
	{
		pp_c<-1
		if(mean(fff_cc[,kk])>mean(fff_cc[,kk+1]))
		{
			pp_c<-t.test(fff_cc[,kk],fff_cc[,kk+1])$p.value
		}
		pp<-c(pp,pp_c)
	}	
	pp_all[[j]]<-pp
	cc<-which(pp<1e-3)
	if(length(cc)>0)
	{
		cc_all<-c(cc_all,cc[length(cc)])	
	}
	if(length(cc)==0)
	{
		cc_all<-c(cc_all,0)	
	}
}
print("BCV done! Compute 1 rank bases and CORS!")
bbb_all<-c()
for(j in 1:length(tg_1_rank_markers))
{
	aaa<-data_CORS_cancer[rownames(tg_1_rank_markers[[j]]),]
	if(cc_all[j]>0)
	{
		bbb<-t(svd(aaa)$v)[1:cc_all[j],]
		if(cc_all[j]==1)
		{
			bbb<-t(as.matrix(t(svd(aaa)$v)[1:cc_all[j],]))
		}
		rownames(bbb)<-c(paste(tg_key_c,j,1:cc_all[j],sep="__"))
		bbb_all<-rbind(bbb_all,bbb)
	}
}
tg_genes_all_c<-c()
for(j in 1:length(tg_1_rank_markers))
{
	tg_genes_all_c<-c(tg_genes_all_c,rownames(tg_1_rank_markers[[j]]))
}
tg_genes_all_c<-sort(unique(tg_genes_all_c))
aaa<-data_CORS_cancer[intersect(rownames(data_CORS_cancer),unique(rownames(immune_cell_uni_table))),]
ccc<-c()
for(i in 1:nrow(bbb_all))
{
	ttt_c<-bbb_all[i,]%*%t(bbb_all[i,])/sum((bbb_all[i,])^2)
	ccc<-cbind(ccc,RMSE_row(aaa%*%ttt_c)/RMSE_row(aaa))
	if(i%%50==1)
	{
		print(c(i,nrow(bbb_all)))
	}
}
colnames(ccc)<-rownames(bbb_all)
print(paste("1 rank bases and CORS done!Compute top",hcutn,"bases"))
tg_F<-paste(tg_key_c,"_bases_correlation.pdf",sep="")
pdf(tg_F)
library(gplots)
colors = c(-100:100)/100
my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n =200)
bbb_all0<-bbb_all[names(which(apply(ccc,2,max)>0.4)),]
heatmap.2(cor(t(bbb_all0)),Rowv=T,Colv =T,scale="none",main=paste(tg_key_c,"\nall selected bases"),
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(10,10),cexRow=0.5,cexCol=0.5)
stat_c<-c()
stat_p_c<-c()
for(j in 1:length(tg_1_rank_markers))
{
	ddd<-apply(tg_1_rank_markers[[j]][,-c(1:4)],2,mean)
	#print(names(which(ddd==max(ddd))))
	#print(j)
	stat_p_c<-c(stat_p_c,names(which(ddd==max(ddd)))[1])
	stat_c<-c(stat_c,max(ddd))
}
names(stat_c)<-stat_p_c	
h<-hclust(dist(bbb_all0))
fff<-cutree(h,hcutn)
fff0<-unique(fff)
tg_trunc_list<-list()
tg_trunc_list_info<-list()
tg_trunc_RS_bases<-c()
for(j in 1:length(fff0))
{	
	tg_ccc<-names(which(fff==fff0[j]))
	tg_trunc_list[[j]]<-tg_ccc
	tg_info_c<-c()
	for(i in 1:length(tg_ccc))
	{
		hhh<-as.numeric(as.character(unlist(strsplit(tg_ccc[i],"__"))[2]))
		tg_info_c<-c(tg_info_c,paste(tg_ccc[i],names(stat_c)[hhh],sep="|"))
	}
	tg_trunc_list_info[[j]]<-tg_info_c
	if(length(tg_ccc)==1)
	{
		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0[tg_ccc,])
	}
	else
	{
		bbb_all0_c<-bbb_all0[tg_ccc,]
		bbb_all0_cc<-t(svd(bbb_all0_c)$v)[1,]
		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0_cc)
	}
}
rownames(tg_trunc_RS_bases)<-1:nrow(tg_trunc_RS_bases)
list_ccc<-list(tg_trunc_RS_bases,tg_trunc_list,tg_trunc_list_info,bbb_all0,stat_c)
names(list_ccc)<-c("RowBase_trunc","RowBase_trunc_ID","RowBase_trunc_info","RowBase_all","RowBase_all_celltype_info")	
heatmap.2(cor(t(tg_trunc_RS_bases)),Rowv=T,Colv =T,scale="none",main=paste(tg_key_c,"\ntruncated bases",hcutn),
col=my_palette,breaks=colors,density.info="none",dendrogram="both",
trace="none",margin=c(10,10),cexRow=0.5,cexCol=0.5)
dev.off()
return(list_ccc)
}



scoring_MR_order<-function(aaa,penalty=-2,lbb=-20)
{
	mm<-max(aaa)
	dd<-1:mm
	cc<-rep(penalty,length(dd))
	cc[aaa]<-1
	ccc<-cumsum(cc)
	if(min(ccc)<=lbb)
	{
		trunc_id<-min(which(ccc<=lbb))
		ccc<-ccc[1:trunc_id]
	}
	tg_id<-which(ccc==max(ccc))
	tg_id<-tg_id[length(tg_id)]
	return(tg_id)
}

calculate_growth_rate2<-function(x0,step=20)
{
        gr_all<-rep(0,length(x0))
        for(i in 1:length(x0))
        {
                sid<-max(1,i-step)
                eid<-i
                tg_id<-c(sid:eid)
                l<-eid-sid+1
                gr_all[i]<-(x0[eid]-x0[sid])/l
        }
        return(gr_all)
}

BCV_screen_top_bases<-function(data_CORS_cancer,tg_key_c=tg_key_c,tg_1_rank_markers=tg_1_rank_markers,immune_cell_uni_table=immune_cell_uni_table0,hcutn=40)
{
library(bcv)
print("start BCV")
pp_all<-list()
cc_all<-c()
for(j in 1:length(tg_1_rank_markers))
{
        aaa<-data_CORS_cancer[c(names(tg_1_rank_markers)[j],rownames(tg_1_rank_markers[[j]])),]
        x<-aaa
        fff_cc<-c()
        for(kk in 1:10)
        {
                cv_result <- cv.svd.gabriel(x, 3, 3, maxrank = 30)
                fff_cc<-rbind(fff_cc,cv_result$msep)
        }
        pp<-c()
        for(kk in 1:(ncol(fff_cc)-1))
        {
                pp_c<-1
                if(mean(fff_cc[,kk])>mean(fff_cc[,kk+1]))
                {
                        pp_c<-t.test(fff_cc[,kk],fff_cc[,kk+1])$p.value
                }
                pp<-c(pp,pp_c)
        }       
        pp_all[[j]]<-pp
        cc<-which(pp<1e-3)
        if(length(cc)>0)
        {
                cc_all<-c(cc_all,cc[length(cc)])        
        }
        if(length(cc)==0)
        {
                cc_all<-c(cc_all,0)     
        }
}
print("BCV done! Compute 1 rank bases and CORS!")
bbb_all<-c()
for(j in 1:length(tg_1_rank_markers))
{
        aaa<-data_CORS_cancer[rownames(tg_1_rank_markers[[j]]),]
        if(cc_all[j]>0)
        {
                bbb<-t(svd(aaa)$v)[1:cc_all[j],]
                if(cc_all[j]==1)
                {
                        bbb<-t(as.matrix(t(svd(aaa)$v)[1:cc_all[j],]))
                }
                rownames(bbb)<-c(paste(tg_key_c,j,1:cc_all[j],sep="__"))
                bbb_all<-rbind(bbb_all,bbb)
        }
}
tg_genes_all_c<-c()
for(j in 1:length(tg_1_rank_markers))
{
        tg_genes_all_c<-c(tg_genes_all_c,rownames(tg_1_rank_markers[[j]]))
}
tg_genes_all_c<-sort(unique(tg_genes_all_c))
aaa<-data_CORS_cancer[intersect(rownames(data_CORS_cancer),unique(rownames(immune_cell_uni_table))),]
ccc<-c()
for(i in 1:nrow(bbb_all))
{
        ttt_c<-bbb_all[i,]%*%t(bbb_all[i,])/sum((bbb_all[i,])^2)
        ccc<-cbind(ccc,RMSE_row(aaa%*%ttt_c)/RMSE_row(aaa))
        if(i%%50==1)
        {
                print(c(i,nrow(bbb_all)))
        }
}
colnames(ccc)<-rownames(bbb_all)
print(paste("1 rank bases and CORS done! Start screen bases"))
stat_c<-c()
stat_p_c<-c()
for(j in 1:length(tg_1_rank_markers))
{
        ddd<-apply(tg_1_rank_markers[[j]][,-c(1:4)],2,mean)
        #print(names(which(ddd==max(ddd))))
        #print(j)
        stat_p_c<-c(stat_p_c,names(which(ddd==max(ddd)))[1])
        stat_c<-c(stat_c,max(ddd))
}
names(stat_c)<-stat_p_c 
bbb_all0<-bbb_all[names(which(apply(ccc,2,max)>0.4)),]
h<-hclust(dist(bbb_all0))
tg_trunc_list<-list()
tg_trunc_list_old<-list()
st<-0
trunc_RS_bases_all<-list()
RMSE_IM<-list()
hclust_info<-list()
for(i in 1:hcutn)
{
	fff<-cutree(h,i)
	hclust_info[[i]]<-fff
	fff0<-unique(fff)
	tg_trunc_list_old<-tg_trunc_list
	tg_trunc_list<-list()
	#tg_trunc_list_info<-list()
	tg_trunc_RS_bases<-c()
	RMSE_table_c<-c()
	if(st==0)
	{
		for(j in 1:length(fff0))
		{
			tg_ccc<-names(which(fff==fff0[j]))
			tg_trunc_list[[j]]<-tg_ccc
			if(length(tg_ccc)==1)
        		{
                		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0[tg_ccc,])
        		}
        		else
        		{
               		bbb_all0_c<-bbb_all0[tg_ccc,]
                		bbb_all0_cc<-t(svd(bbb_all0_c)$v)[1,]
				aaa0<-cor(t(t(bbb_all0_cc)),t(bbb_all0_c))
				aaa0<-aaa0[which(abs(aaa0)>0.6)]
				if(mean(aaa0)<0)
				{
					bbb_all0_cc<--bbb_all0_cc
				}
               		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0_cc)
        		}
		}
		rownames(tg_trunc_RS_bases)<-1:nrow(tg_trunc_RS_bases)
		trunc_RS_bases_all[[i]]<-tg_trunc_RS_bases
		ccc<-c()
		for(j in 1:nrow(tg_trunc_RS_bases))
		{
        		ttt_c<-tg_trunc_RS_bases[j,]%*%t(tg_trunc_RS_bases[j,])/sum((tg_trunc_RS_bases[j,])^2)
        		ccc<-cbind(ccc,RMSE_row(aaa%*%ttt_c)/RMSE_row(aaa))
        	}
		RMSE_IM[[i]]<-ccc
		st<-1
	}
	else
	{
		ccc_old<-ccc
		for(j in 1:length(fff0))
		{
			tg_ccc<-names(which(fff==fff0[j]))
			tg_trunc_list[[j]]<-tg_ccc
			if(length(tg_ccc)==1)
        		{
                		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0[tg_ccc,])
        		}
        		else
        		{
               		bbb_all0_c<-bbb_all0[tg_ccc,]
                		bbb_all0_cc<-t(svd(bbb_all0_c)$v)[1,]
				aaa0<-cor(t(t(bbb_all0_cc)),t(bbb_all0_c))
				aaa0<-aaa0[which(abs(aaa0)>0.6)]
				if(mean(aaa0)<0)
				{
					bbb_all0_cc<--bbb_all0_cc
				}
               		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0_cc)
              	}
		}
		rownames(tg_trunc_RS_bases)<-1:nrow(tg_trunc_RS_bases)
		trunc_RS_bases_all[[i]]<-tg_trunc_RS_bases
		perturb_id<-hcut_update(tg_trunc_list,tg_trunc_list_old)
		ccc<-c()
		for(j in 1:nrow(tg_trunc_RS_bases))
		{
			if(sum(j==perturb_id[[2]])==0)
			{
				ccc<-cbind(ccc,ccc_old[,perturb_id[[3]][j]])
			}
			else
			{
				ttt_c<-tg_trunc_RS_bases[j,]%*%t(tg_trunc_RS_bases[j,])/sum((tg_trunc_RS_bases[j,])^2)
				ccc<-cbind(ccc,RMSE_row(aaa%*%ttt_c)/RMSE_row(aaa))
				print(c(i,j))


			}
        	}
		RMSE_IM[[i]]<-ccc
	}
}
N<-hcutn
tg_ids<-seq(hcutn+1,nrow(bbb_all0),by=30)
for(i in 1:length(tg_ids))
{
	N<-N+1
	print(c(N,tg_ids[i]))
	fff<-cutree(h,tg_ids[i])
	hclust_info[[N]]<-fff
	fff0<-unique(fff)
	tg_trunc_list<-list()
	tg_trunc_RS_bases<-c()
	RMSE_table_c<-c()
	st<-0
	if(st==0)
	{
		for(j in 1:length(fff0))
		{
			tg_ccc<-names(which(fff==fff0[j]))
			tg_trunc_list[[j]]<-tg_ccc
			if(length(tg_ccc)==1)
        		{
                		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0[tg_ccc,])
        		}
        		else
        		{
               		bbb_all0_c<-bbb_all0[tg_ccc,]
                		bbb_all0_cc<-t(svd(bbb_all0_c)$v)[1,]
				aaa0<-cor(t(t(bbb_all0_cc)),t(bbb_all0_c))
				aaa0<-aaa0[which(abs(aaa0)>0.6)]
				if(mean(aaa0)<0)
				{
					bbb_all0_cc<--bbb_all0_cc
				}
               		tg_trunc_RS_bases<-rbind(tg_trunc_RS_bases,bbb_all0_cc)
        		}
		}
		rownames(tg_trunc_RS_bases)<-1:nrow(tg_trunc_RS_bases)
		trunc_RS_bases_all[[N]]<-tg_trunc_RS_bases
		ccc<-c()
		for(j in 1:nrow(tg_trunc_RS_bases))
		{
        		ttt_c<-tg_trunc_RS_bases[j,]%*%t(tg_trunc_RS_bases[j,])/sum((tg_trunc_RS_bases[j,])^2)
        		ccc<-cbind(ccc,RMSE_row(aaa%*%ttt_c)/RMSE_row(aaa))
        	}
		RMSE_IM[[N]]<-ccc
	}
}
ccc<-list(trunc_RS_bases_all,RMSE_IM,hclust_info)
names(ccc)<-c("Bases_all","RMSE_all","clust_info_all")
return(ccc)
}

hcut_update<-function(tg_trunc_list,tg_trunc_list_old)
{
	out_list<-rep(0,length(tg_trunc_list))
	for(i in 1:length(tg_trunc_list))
	{
		for(j in 1:length(tg_trunc_list_old))
		{
			if(length(intersect(tg_trunc_list[[i]],tg_trunc_list_old[[j]]))>0)
			{
				out_list[i]<-j
			}
		}
	}
	tg_old_id<-as.numeric(names(which(table(out_list)==2)))
	tg_new_id<-which(out_list==tg_old_id)
	return(list(tg_old_id,tg_new_id,out_list))
}

plot_simple2<-function(MR_list,tg_gene="",ll=200)
{
          bbb<-MR_list[,3]
	    ll0<-min(ll,nrow(MR_list))
          ccc<-MR_list[,-c(1:4)]
	    ccc1<-apply(ccc,2,mean)
	    tg_id<-which(ccc1==max(ccc1))[1]
 	    plot(bbb,type="b",pch=16,xlim=c(1,ll0),ylim=c(0,2),lwd=4,col="lightblue",main=paste(tg_gene,"\n",colnames(ccc)[tg_id],sep=""))
          points((cumsum(ccc[,tg_id])/c(1:nrow(ccc))),pch=16,type="b",lwd=4,col="palegreen")
          bbb2<-(cumsum(ccc[,tg_id])/c(1:nrow(ccc)))
 
          bbb3<-MR_list[,2]
          points(bbb3,pch=16,type="b",lwd=4,col="violet",ylim=c(-1,2))

          bbb4<-MR_list[,4]
          points(bbb4,pch=16,type="b",lwd=4,col="orange")
	    abline(h=0.6,col="blue",lty=2)
	    abline(h=0.7,col="red",lty=2)
	    abline(h=0.8,col="blue",lty=2)
	    abline(h=1,col="black",lty=2)
	    legend(ll0*0.75,1.8,c("MR_sorted","Cell_Enrich","MR_growth","Corr"),lwd=4,col=c("lightblue","palegreen","violet","orange"))
}

MAP_GPL570_genes<-function(tg_cluster_genes)
{
	tg_cluster_genes0<-tg_cluster_genes
	for(i in 1:length(tg_cluster_genes))
	{
		tg_cluster_genes0[[i]]<-GPL570_id_symbol0[tg_cluster_genes[[i]],2]
	}
	return(tg_cluster_genes0)
}


cut_modules<-function(tg_list,cutn=20)
{
	tg_list0<-tg_list
	for(i in 1:length(tg_list0))
	{
		ccc<-tg_list[[i]]
		nn<-min(length(ccc),cutn)
		tg_list0[[i]]<-ccc[1:nn]
	}
	names(tg_list0)<-names(tg_list)
	return(tg_list0)
}

MAP_GPL570_genes2<-function(tg_cluster_genes)
{
	tg_cluster_genes0<-tg_cluster_genes
	for(i in 1:length(tg_cluster_genes))
	{
		tg_cluster_genes0[[i]]<-sort(unique(GPL570_id_symbol0[tg_cluster_genes[[i]],2]))
	}
	return(tg_cluster_genes0)
}



R1_list_filtering_step1<-function(list_c2,data_CORS_cancer,data01,data0,cutn0=15,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
tg_1_rank_markers0<-list_c2[[1]]
tg_m_names0<-list_c2[[2]]
RMSE_on_CORS_cancer_c<-RMSE_row(data_CORS_cancer)
tg_marker_lists<-list()
for(i in 1:length(list_c2[[1]]))
{
	ccc<-c(names(list_c2[[1]])[i],rownames(list_c2[[1]][[i]]))
	tg_marker_lists[[i]]<-ccc
}
names(tg_marker_lists)<-names(list_c2[[1]])

tg_marker_lists<-cut_modules(tg_marker_lists,cutn=cutn0)
names(tg_marker_lists)<-names(list_c2[[1]])
tg_genes_all<-c()
for(i in 1:length(tg_marker_lists))
{
	tg_genes_all<-c(tg_genes_all,tg_marker_lists[[i]])
}
tg_genes_all<-unique(sort(tg_genes_all))

data_c1<-data_CORS_cancer[tg_genes_all,]
data_c2<-data01[tg_genes_all,]
data_c3<-data0[tg_genes_all,]
cluster_n<-c()
for(i in 1:length(tg_marker_lists))
{
	cluster_n<-c(cluster_n,length(tg_marker_lists[[i]]))
}
RMSE_CORSC_hub<-c()
for(i in 1:length(tg_marker_lists))
{
	RMSE_CORSC_hub<-c(RMSE_CORSC_hub,RMSE_on_CORS_cancer_c[names(tg_marker_lists)[i]])
}
RMSE_CORSC_all<-c()
for(i in 1:length(tg_marker_lists))
{
	RMSE_CORSC_all<-c(RMSE_CORSC_all,mean(RMSE_on_CORS_cancer_c[tg_marker_lists[[i]]]))
}
Top_cell_proportion<-c()
stat_p_c<-c()
for(j in 1:length(tg_1_rank_markers0))
{
        ddd<-apply(tg_1_rank_markers0[[j]][,-c(1:4)],2,mean)
        stat_p_c<-c(stat_p_c,names(which(ddd==max(ddd)))[1])
        Top_cell_proportion<-c(Top_cell_proportion,max(ddd))
}
names(Top_cell_proportion)<-stat_p_c
cluster_stat<-cbind(cluster_n,Top_cell_proportion,RMSE_CORSC_hub,RMSE_CORSC_all)

stat_CORS_cancer<-evaluate_rank1_markers(data_c1,tg_marker_lists)
stat_FPKM_lognew_normalized<-evaluate_rank1_markers(data_c2,tg_marker_lists)
stat_FPKM_lognew<-evaluate_rank1_markers(data_c3,tg_marker_lists)

R1_markers_f1<-tg_marker_lists
R1_stat_all_c<-list(cluster_stat,Top_cell_proportion,R1_markers_f1,stat_CORS_cancer,stat_FPKM_lognew_normalized,stat_FPKM_lognew)
print("filter 1 and stat done!")

pp_all<-c()
for(i in 1:length(R1_markers_f1))
{
	pp<-sum(BCV_ttest(data_CORS_cancer[R1_markers_f1[[i]],],maxrank0=20)<0.001)	
	pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f1<-pp_all

#ccc<-MAP_GPL570_genes2(R1_markers_f1)
#names(ccc)<-names(Top_cell_proportion)
#table(names(Top_cell_proportion))

cut1<-cut10
ccc<-compute_min_jaccard(R1_markers_f1)
ccc0<-ccc>cut1

stat_cc<-c(1:nrow(ccc0))
names(stat_cc)<-1:nrow(ccc0)
for(i in 1:nrow(ccc0))
{
	for(j in 1:ncol(ccc0))
	{
		if((i<j)&(ccc0[i,j]>0))
		{
			nn<-max(i,j)
			stat_cc[which(stat_cc==i)]<-nn
			stat_cc[which(stat_cc==j)]<-nn
		}
	}
}
table(stat_cc)
tg_ccc<-unique(stat_cc)
R1_marker_list_f2<-list()
N<-0
for(i in 1:length(tg_ccc))
{
	N<-N+1	
	tg_ids<-as.numeric(names(stat_cc)[which(stat_cc==tg_ccc[i])])
	ccc<-c()
	for(j in 1:length(tg_ids))
	{
		ccc<-c(ccc,R1_markers_f1[[tg_ids[j]]])
	}
	ccc<-unique(ccc)
	R1_marker_list_f2[[N]]<-ccc
}
print("filter 2 done!")
#R1_marker_list_f2<-filter_1rank_genes(R1_marker_list_f2,data_CORS_cancer,cut0=0.7)
ccc<-c()
nn<-c()
for(i in 1:length(R1_marker_list_f2))
{	
	ddd<-apply(immune_cell_uni_table[R1_marker_list_f2[[i]],],2,mean)
	ccc<-rbind(ccc,ddd)
	nn<-c(nn,colnames(immune_cell_uni_table)[which(ddd==max(ddd))[1]])
}
#eee<-MAP_GPL570_genes2(R1_marker_list_f2)
names(R1_marker_list_f2)<-nn
#names(eee)<-nn
pp_all<-c()
for(i in 1:length(R1_marker_list_f2))
{
	pp<-sum(BCV_ttest(data_CORS_cancer[R1_marker_list_f2[[i]],],maxrank0=20)<0.001)	
	pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f2<-pp_all

Filter_1_result_list<-list(R1_markers_f1,R1_stat_all_c,pp_R1_marker_list_f1,
R1_marker_list_f2,pp_R1_marker_list_f2)
names(Filter_1_result_list)<-c("R1_markers_f1","R1_stat_all_c","pp_R1_marker_list_f1",
"R1_marker_list_f2","pp_R1_marker_list_f2")
return(Filter_1_result_list)
}

R1_list_filtering_step2_BCV<-function(R1_marker_list_f2,data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,tg_top_markers=tg_top_markers_GPL570)
{

cut1<-cut10
#library(gplots)
#colors = c(0:100)/100
#my_palette <- grDevices::colorRampPalette(c("white", "blue"))(n =100)

ccc<-compute_min_jaccard(R1_marker_list_f2)
ccc0<-ccc>cut1
names(R1_marker_list_f2)<-1:length(R1_marker_list_f2)
print("filter 3 BCV start!")
module_merge_info<-compute_BCV_1fold_test(R1_marker_list_f2,data_CORS_cancer,ccc0)
print("filter 3 BCV done!")

tg_cluster_genes<-module_merge_info[[2]]
pp_all<-c()
for(i in 1:length(tg_cluster_genes))
{
	pp<-sum(BCV_ttest(data_CORS_cancer[tg_cluster_genes[[i]],],maxrank0=20)<0.001)	
	pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-pp_all
stat0_all<-c()
for(i in 1:length(tg_cluster_genes))
{
	stat0<-c()
	for(j in 1:length(tg_top_markers))
	{
		#stat0<-c(stat0,length(intersect(tg_top_markers[[j]],tg_cluster_genes[[i]]))/length(tg_cluster_genes[[i]]))
		#stat0<-c(stat0,length(intersect(tg_top_markers[[j]],tg_cluster_genes[[i]]))/length(tg_cluster_genes[[i]]))
	}
	stat0<-apply(sqrt(immune_cell_uni_table[tg_cluster_genes[[i]],]),2,mean)
	stat0_all<-rbind(stat0_all,stat0)
}
colnames(stat0_all)<-names(tg_top_markers)
rownames(stat0_all)<-1:nrow(stat0_all)
pp<-apply(stat0_all,1,max)
pn<-c()
for(i in 1:nrow(stat0_all))
{
	pn<-c(pn,colnames(stat0_all)[which(stat0_all[i,]==max(stat0_all[i,]))[1]])
}
names(pp)<-pn
names(tg_cluster_genes)<-pn
#MAP_GPL570_genes2(tg_cluster_genes)


	tg_cluster_bases<-compute_base2(tg_cluster_genes,data_CORS_cancer)
	rownames(tg_cluster_bases)<-pn
	cor0<-cor(t(tg_cluster_bases))
	rownames(cor0)<-pn
	colnames(cor0)<-pn
	cor1<-compute_adjusted_cor_higher(tg_cluster_genes,data_CORS_cancer,qant=0.5)
	rownames(cor1)<-pn
	colnames(cor1)<-pn
	cor2<-compute_adjusted_cor_lower(tg_cluster_genes,data_CORS_cancer,qant=0.5)
	rownames(cor2)<-pn
	colnames(cor2)<-pn
	cor3<-compute_adjusted_cor_all(tg_cluster_genes,data_CORS_cancer,qanth=0.5,qantl=0.5)
	rownames(cor3)<-pn
	colnames(cor3)<-pn
	#library(gplots)
	#colors = c(-100:100)/100
	#my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n =200)

	#heatmap.2(cor0,Rowv=T,Colv =T,scale="none",main="original",
	#col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	#trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

	#h<-heatmap.2(cor1,Rowv=T,Colv =T,scale="none",main="high",
	#col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	#trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

	#heatmap.2(cor3,Rowv=T,Colv =T,scale="none",main="low",
	#col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	#trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

	#heatmap.2(cor3^2,Rowv=T,Colv =T,scale="none",main="all",
	#col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	#trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

	Filter_2BCV_result_list<-list(tg_cluster_genes,pp_R1_marker_list_f3,cor0,cor1,cor2,cor3)
	names(Filter_2BCV_result_list)<-c("tg_cluster_genes","pp_R1_marker_list_f3","cor0","cor1","cor2","cor3")
	return(Filter_2BCV_result_list)
}


clean_rank1_module<-function(data_c,module_info,module_rank,st0=8)
{
	N<-0
	nc<-c()
	module_new<-list()
	for(i in 1:length(module_info))
	{
		if(module_rank[i]==1)
		{
			N<-N+1
			nc<-c(nc,names(module_info)[i])
			module_new[[N]]<-module_info[[i]]		
		}
		if(module_rank[i]>1)
		{
			ccc<-module_info[[i]]	
			st<-st0
			rr<-1
			while((rr==1)&(st<=length(ccc)))
			{
				tg_genes<-c(ccc[1:st])
				pp<-BCV_ttest(data_c[tg_genes,],maxrank0=5)
				rr<-sum(pp<0.001)
				st<-st+1
			}
			tg_genes<-tg_genes[-length(tg_genes)]
			if(length(tg_genes)>st0)
			{
				N<-N+1
				nc<-c(nc,names(module_info)[i])
				module_new[[N]]<-tg_genes	
			}
		}
	}
	names(module_new)<-nc
	return(module_new)
}



R1_list_filtering_step3_buildNMF<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
data_CORS_cancer<-data_c
module_info<-R1_filter_step1_results[[4]]
module_rank<-R1_filter_step1_results[[5]]
R1_marker_list_f2.5<-clean_rank1_module(data_c,module_info,module_rank,st0=4)
cut1<-cut10
ccc<-compute_min_jaccard(R1_marker_list_f2.5)
ccc0<-ccc>cut1
names(R1_marker_list_f2.5)<-1:length(R1_marker_list_f2.5)
print("Run 1st round BCV and Jaccard distance for merging bases!")
module_merge_info<-compute_BCV_1fold_test(R1_marker_list_f2.5,data_c,ccc0)
tg_cluster_genes<-module_merge_info[[2]]
pp_all<-c()
for(i in 1:length(tg_cluster_genes))
{
        pp<-sum(BCV_ttest(data_c[tg_cluster_genes[[i]],],maxrank0=20)<0.001)  
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-pp_all
R1_marker_list_f3.5<-clean_rank1_module(data_c,tg_cluster_genes,pp_R1_marker_list_f3,st0=4)
cut1=0.8
ccc<-compute_min_jaccard(R1_marker_list_f3.5)
ccc0<-ccc>cut1
names(R1_marker_list_f3.5)<-1:length(R1_marker_list_f3.5)
print("Constructing Linking Graph!")
BCV_new_stat<-BCV_1fold_test_fast(marker_list=R1_marker_list_f3.5,data_c=data_c)
ccc<-sign((BCV_new_stat[[1]]==1)+(BCV_new_stat[[3]]==1)+(BCV_new_stat[[2]]==1)+ccc0)
diag(ccc)<-1
heatmap.2(ccc
,Rowv=T,Colv =T,scale="none",main="",
        col=my_palette,breaks=colors,density.info="none",dendrogram="both",
        trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)
graph_links<-ccc
marker_list<-R1_marker_list_f3.5
g1 <- graph_from_adjacency_matrix(ccc)
largest_cliques(g1)
cg<-clusters(g1)
tg_ccc<-as.numeric(names(which(cg$membership==5)))
tg_clusters_all<-unique(cg$membership)
cluster_number<-c()
for(i in 1:length(tg_clusters_all))
{
	cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
}
R1_marker_merged_gc<-list()
for(i in 1:length(tg_clusters_all))
{
	ccc<-c()
	tgs<-which(cg$membership==tg_clusters_all[i])
	for(j in 1:length(tgs))
	{
		ccc<-c(ccc,marker_list[[tgs[j]]])
	}
	R1_marker_merged_gc[[i]]<-unique(ccc)
}
tg1c<-which(cluster_number==1)
tgoc<-which(cluster_number>1)
oc_1c_cor<-matrix(0,length(tg1c),length(tgoc))
rownames(oc_1c_cor)<-tg1c
colnames(oc_1c_cor)<-tgoc
for(i in 1:length(tg1c))
{
	for(j in 1:length(tgoc))
	{
		tg1<-R1_marker_merged_gc[[tg1c[i]]]
		tg2<-R1_marker_merged_gc[[tgoc[j]]]
		k1<-sum(BCV_ttest(data_c[tg2,],rounds=40)<0.001)
		k2<-sum(BCV_ttest(data_c[tg1,],rounds=40)<0.001)
		k3<-sum(BCV_ttest(data_c[union(tg1,tg2),],rounds=40)<0.001)
		cc1<-apply(immune_cell_uni_table[tg2,],2,mean)
		cc2<-apply(immune_cell_uni_table[tg1,],2,mean)
		dd1<-c(cc2[which(cc2==max(cc2))[1]]-cc1[which(cc2==max(cc2))[1]])
		dd2<-c(cc1[which(cc1==max(cc1))[1]]-cc2[which(cc1==max(cc1))[1]])
		#print(c(i,j))
		#print(c(dd1,dd2))
		if((k3<=k1)&(dd1<0.5)&(dd2<0.5))
		{
			oc_1c_cor[i,j]<-1
		}
	}
}
tg1c_new<-c()
if(sum(apply(oc_1c_cor,1,sum)==0)>0)
{
	tg1c_new<-tg1c[which(apply(oc_1c_cor,1,sum)==0)]
}
tg1c_merge<-setdiff(tg1c,tg1c_new)
All_genes<-c()
for(i in 1:length(marker_list))
{
	All_genes<-c(All_genes,marker_list[[i]])
}
All_genes<-sort(unique(All_genes))
NMF_indi_all<-c()
print("Working on Linking graph for NMF constraints!")
for(ii in 1:length(tgoc))
{
	print(c(ii,length(tgoc)))
	tg1c_m_c<-tg1c[which(oc_1c_cor[,ii]==1)]
	tgs<-which(cg$membership==tgoc[ii])
	if(length(tg1c_m_c)>0)
	{
		for(j in 1:length(tg1c_m_c))
		{
			tgs<-c(tgs,which(cg$membership==tg1c_m_c[j]))
		}
	}
	tg_M_all<-c()
	for(j in 1:length(tgs))
	{
		tg_M_all<-c(tg_M_all,marker_list[[tgs[j]]])
	}
	tg_M_all<-unique(tg_M_all)
	pp<-c()
	for(k in 1:21)
	{
		pp<-c(pp,sum(BCV_ttest(data_c[tg_M_all,],rounds=40,maxrank0=20)<0.01))
	}
	pp0<-floor(median(pp))
	Base_all<-c()
	for(j in 1:length(tgs))
	{
		Base_all<-rbind(Base_all,apply(data_c[marker_list[[tgs[j]]],],2,mean))
	}
	NMF_label_c<-c()
	if(pp0<=1)
	{
		NMF_label_c<-rep(1,length(tg_M_all))
		names(NMF_label_c)<-tg_M_all
		NMF_label_c<-as.matrix(NMF_label_c)
		colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste(tgoc[ii],"All1",sep="_")
	}
	if(pp0>1)
	{
	stat_all<-rep(1,length(tgs))
	names(stat_all)<-tgs
	if(length(tgs)>=3)
	{
	ddd<-combn(length(tgs),3)
	comp_ids<-c()
	stat_all<-rep(0,length(tgs))
	names(stat_all)<-tgs
	for(i in 1:ncol(ddd))
	{
		tg_genes<-unique(c(marker_list[[tgs[ddd[1,i]]]],marker_list[[tgs[ddd[2,i]]]],marker_list[[tgs[ddd[3,i]]]]))
		pp<-sum(BCV_ttest(data_c[tg_genes,],rounds=20,maxrank0=5)<0.01)
		if(pp==2)
		{
			comp_ids<-c(comp_ids,i)
			#print(i)
			lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
			if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==1))
			{
				stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
				stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
			}
			if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[2]==1))
			{
				stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
				stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
			}
			if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==-1))
			{
				stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
				stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
			}
		}
	}
	}
	ccc<-stat_all[order(-stat_all)]
	ccc<-ccc[which(ccc!=0)]
	specific_cell_marker_ids<-c()
	if(length(ccc)>0)
	{
	specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[1]))
	for(j in 1:length(ccc))
	{
		if(sum(graph_links[specific_cell_marker_ids,as.numeric(names(ccc)[j])])==0)
		{
			specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[j]))
		}
	}
	if(length(specific_cell_marker_ids)>pp0)
	{
		specific_cell_marker_ids<-specific_cell_marker_ids[1:pp0]
	}
	NMF_label_c<-c()
	if(length(specific_cell_marker_ids)>1)
	{
		NMF_label_c<-matrix(0,length(tg_M_all),length(specific_cell_marker_ids))
		rownames(NMF_label_c)<-tg_M_all
		colnames(NMF_label_c)<-specific_cell_marker_ids
		tg_specific_gene<-c()
		for(i in 1:length(specific_cell_marker_ids))
		{
			tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
			NMF_label_c[marker_list[[specific_cell_marker_ids[i]]],i]<-1
		}
		tg_all_genes<-setdiff(tg_M_all,unique(tg_specific_gene))
		NMF_label_c[tg_all_genes,]<-NMF_label_c[tg_all_genes,]+1
	}
	if(length(specific_cell_marker_ids)==1)
	{
		NMF_label_c<-matrix(0,length(tg_M_all),length(specific_cell_marker_ids))
		rownames(NMF_label_c)<-tg_M_all
		colnames(NMF_label_c)<-specific_cell_marker_ids
		tg_specific_gene<-c()
		for(i in 1:length(specific_cell_marker_ids))
		{
			tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
			NMF_label_c[marker_list[[specific_cell_marker_ids[i]]],i]<-1
		}
		NMF_label_c<-cbind(NMF_label_c,rep(1,length(tg_M_all)))
		colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste(tgoc[ii],"All1",sep="_")
	}
	}
	if((length(specific_cell_marker_ids)==0)&pp<3)
	{
		NMF_label_c<-rep(1,length(tg_M_all))
		names(NMF_label_c)<-tg_M_all
		NMF_label_c<-as.matrix(NMF_label_c)
		colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1",sep="_")
	}
	if((length(specific_cell_marker_ids)==0)&pp>=3)
	{
		NMF_label_c<-c()
		NMF_label_c0<-rep(0,length(tg_M_all))
		names(NMF_label_c0)<-tg_M_all
		for(j in 1:length(tgs))
		{
			NMF_label_ccc<-NMF_label_c0
			NMF_label_ccc[marker_list[[tgs[j]]]]<-1
			NMF_label_c<-cbind(NMF_label_c,NMF_label_ccc)
		}
		colnames(NMF_label_c)<-paste("tgoc",tgs,"ss1",sep="_")
	}
	}
	NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
	colnames(NMF_label_c_all)<-colnames(NMF_label_c)
	rownames(NMF_label_c_all)<-All_genes
	NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
	NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
}
if(length(tg1c_new)>0)
{
for(ii in 1:length(tg1c_new))
{
	tgs<-which(cg$membership==tg1c_new[ii])
	tg_M_all<-marker_list[[tgs]]
	NMF_label_c<-rep(1,length(tg_M_all))
	names(NMF_label_c)<-tg_M_all
	NMF_label_c<-as.matrix(NMF_label_c)
	colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tg1c",tg1c_new[ii],"All1",sep="_")
	NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
	colnames(NMF_label_c_all)<-colnames(NMF_label_c)
	rownames(NMF_label_c_all)<-All_genes
	NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
	NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
}
}
R1_filter_step3_results<-list(NMF_indi_all,R1_marker_list_f2.5,R1_marker_list_f3.5,graph_links,g1,cg)
names(R1_filter_step3_results)<-c("NMF_indi_all","R1_marker_list_f2.5","R1_marker_list_f3.5","graph_links","constructed graph","graph cluster")
return(R1_filter_step3_results)
}


BCV_1fold_test_fast<-function(marker_list,data_c,pp=0.001)
{
stat<-matrix(0,length(marker_list),length(marker_list))
rownames(stat)<-names(marker_list)
colnames(stat)<-names(marker_list)
stat_i<-stat
stat_j<-stat
stat_k<-stat
for(i in 1:length(marker_list))
{
                print(c(i,length(marker_list)))
                for(j in 1:length(marker_list))
                {
                        if(i<j)
                        {       
                                ccc1<-data_c[marker_list[[i]],]
                                ccc2<-data_c[marker_list[[j]],]
                                ccc<-rbind(ccc1,ccc2)
                                ccc<-ccc[unique(rownames(ccc)),]
                                tg_samples1<-which(apply(ccc1,2,mean)>quantile(apply(ccc1,2,mean),0.5))
                                tg_samples2<-which(apply(ccc2,2,mean)>quantile(apply(ccc2,2,mean),0.5))
                                i12<-sum(BCV_ttest(ccc,maxrank0=10)<pp)
                                j12<-sum(BCV_ttest(ccc[,tg_samples1],maxrank0=10)<pp)
                                k12<-sum(BCV_ttest(ccc[,tg_samples2],maxrank0=10)<pp)
                                stat_i[i,j]<-stat_i[j,i]<-i12
                                stat_j[i,j]<-stat_j[j,i]<-j12
                                stat_k[i,j]<-stat_k[j,i]<-k12
                                }
                    }
}
return(list( stat_i, stat_j, stat_k))
}


compute_CompRowspace_NN<-function(data_cc=tg_data,module_RB=module_selected[[4]],ROUNDS=3)
{
	data_CORS_cancer<-data_cc
	for(ii in 1:ROUNDS)
	{
	ccc<-cor(t(data_CORS_cancer),t(module_RB))
	if(nrow(module_RB)>1)
	{
		ddd<-apply(ccc,1,order)
		eee<-ddd[nrow(ddd),]
		eee<-eee*(apply(ccc,1,max)>0)
		tg_ids<-sort(setdiff(unique(eee),0))
		if(length(tg_ids)>0)
		{
			for(i in 1:length(tg_ids))
			{
				ttt_ccc<-module_RB[tg_ids[i],]%*%t(module_RB[tg_ids[i],])/sum((module_RB[tg_ids[i],])^2)
				data_CORS_cancer[names(which(eee==tg_ids[i])),]<-data_CORS_cancer[names(which(eee==tg_ids[i])),]-data_CORS_cancer[names(which(eee==tg_ids[i])),]%*%ttt_ccc
			}
		}
	}
	else
	{
		ddd<-apply(ccc,1,order)
		eee<-ddd
		eee<-eee*(as.vector(ccc)>0)
		tg_ids<-sort(setdiff(unique(eee),0))
		if(length(tg_ids)>0)
		{
			for(i in 1:length(tg_ids))
			{
				ttt_ccc<-module_RB[tg_ids[i],]%*%t(module_RB[tg_ids[i],])/sum((module_RB[tg_ids[i],])^2)
				data_CORS_cancer[names(which(eee==tg_ids[i])),]<-data_CORS_cancer[names(which(eee==tg_ids[i])),]-data_CORS_cancer[names(which(eee==tg_ids[i])),]%*%ttt_ccc
			}
		}
	}
	}
	return(data_CORS_cancer)
}


Select_cancer_module<-function(cancer_module,data_cc=data_t,kk=1)
{
	tg_genes_all=rownames(data_cc)
	#load("Cancer_Pathway_selected.RData")
	cc<-c()
	Module_PE<-list()
	for(i in 1:length(cancer_module[[1]]))
	{
		#print(i)
		ccc<-c()
		if((nrow(cancer_module[[1]][[i]])>100))
		{
			ccc<-cor(t(cancer_module[[1]][[i]]))
			tg_ll<-which(apply(ccc<0,1,sum)<(0.3*nrow(ccc)))
			cancer_module[[1]][[i]]<-cancer_module[[1]][[i]][tg_ll,]
			ccc<-calculation_mutation_TCGA_2017(rownames(cancer_module[[1]][[i]])[tg_ll],pathway_list=Cancer_Pathway_selected,all_genes=tg_genes_all)
			Module_PE[[i]]<-ccc
		}
		cc<-c(cc,sum(ccc[,3]<0.05))
		#print(rownames(ccc)[which(ccc[,3]<0.05)])
	}
	dd<-unique(c(order(-cc)))
	N<-0
	tg_module_selected<-list()
	tg_module_PE_selected<-list()
	tg_ids<-c()
	for(i in 1:length(dd))
	{
		if(cc[dd[i]]>0)
		{
			if(N<kk)
			{
				N<-N+1
				tg_module_selected[[N]]<-rownames(cancer_module[[1]][[dd[i]]])
				tg_module_PE_selected[[N]]<-Module_PE[[dd[i]]]
				tg_ids<-c(tg_ids,dd[i])
			}
		}
	}
	module_base1<-c()
	for(i in 1:length(tg_module_selected))
	{
		tg_data_c<-data_cc[intersect(tg_module_selected[[i]],rownames(data_cc)),]
		cc<-svd(tg_data_c)$v[,1]
		ccc<-cor(cc,t(tg_data_c))
		if(mean(ccc)<0)
		{
			cc<--cc
		}
		module_base1<-rbind(module_base1,cc)
	}
	colnames(module_base1)<-colnames(data_cc)
	rownames(module_base1)<-1:nrow(module_base1)
	nn<-min(kk,nrow(module_base1))
	nn<-max(nn,1)
	module_base2<-(t(svd(module_base1)$v)[1:nn,])
	if(nn==1)
	{
		module_base2<-t(as.matrix(module_base2))
	}
	colnames(module_base2)<-colnames(data_cc)
	rownames(module_base2)<-1:nrow(module_base2)
	tg_genes_all<-c()
	for(i in 1:length(tg_module_selected))
	{
		tg_genes_all<-c(tg_genes_all,tg_module_selected[[i]])
	}
	tg_genes_all<-unique(tg_genes_all)
	tg_data_c<-data_cc[intersect(tg_genes_all,rownames(data_cc)),]
	nn<-min(kk,nrow(module_base1))
	nn<-max(nn,1)
	module_base3<-(t(svd(tg_data_c)$v)[1:nn,])
	if(nn==1)
	{
		module_base3<-t(as.matrix(module_base3))
	}
	colnames(module_base3)<-colnames(data_cc)
	rownames(module_base3)<-1:nrow(module_base3)
	module_selected<-list(tg_module_selected,tg_module_PE_selected,tg_ids,module_base1,module_base2,module_base3)
	names(module_selected)<-c("Coexp_module_selected_genes","Coexp_module_selected_PES","Coexp_module_selected_IDs","module_base1","module_base2","module_base3")
	return(module_selected)
}

BCV_1fold_test_fast_diag<-function(marker_list,data_c)
{
stat<-matrix(0,length(marker_list),length(marker_list))
rownames(stat)<-names(marker_list)
colnames(stat)<-names(marker_list)
stat_i<-stat
stat_j<-stat
stat_k<-stat
for(i in 1:length(marker_list))
{
	ccc<-data_c[marker_list[[i]],]
      ccc<-ccc[unique(rownames(ccc)),]
      tg_samples1<-which(apply(ccc,2,mean)>quantile(apply(ccc,2,mean),0.5))
      tg_samples2<-which(apply(ccc,2,mean)>quantile(apply(ccc,2,mean),0.5))
      i12<-sum(BCV_ttest(ccc,maxrank0=10)<0.001)
      j12<-sum(BCV_ttest(ccc[,tg_samples1],maxrank0=10)<0.001)
      k12<-sum(BCV_ttest(ccc[,tg_samples2],maxrank0=10)<0.001)
      stat_i[i,i]<-i12
    	stat_j[i,i]<-j12
     	stat_k[i,i]<-k12
}
return(list( stat_i, stat_j, stat_k))
}

R1_list_filtering_step25_build_LG<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3,BCVpp=0.001)
{
data_CORS_cancer<-data_c
module_info<-R1_filter_step1_results[[4]]
module_rank<-R1_filter_step1_results[[5]]
R1_marker_list_f2.5<-clean_rank1_module(data_c,module_info,module_rank,st0=6)
cut1<-cut10
ccc<-compute_min_jaccard(R1_marker_list_f2.5)
ccc0<-ccc>cut1
names(R1_marker_list_f2.5)<-1:length(R1_marker_list_f2.5)
print("Run 1st round BCV and Jaccard distance for merging bases!")
module_merge_info<-compute_BCV_1fold_test(R1_marker_list_f2.5,data_c,ccc0)
print("First BCV done!")
tg_cluster_genes<-module_merge_info[[2]]
pp_all<-c()
for(i in 1:length(tg_cluster_genes))
{
        pp<-sum(BCV_ttest(data_c[tg_cluster_genes[[i]],],maxrank0=20)<0.001)  
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-pp_all
R1_marker_list_f3.5<-clean_rank1_module(data_c,tg_cluster_genes,pp_R1_marker_list_f3,st0=6)
cut1=cut10
ccc<-compute_min_jaccard(R1_marker_list_f3.5)
ccc0<-ccc>cut1
rownames(ccc0)<-1:nrow(ccc0)
colnames(ccc0)<-1:ncol(ccc0)

names(R1_marker_list_f3.5)<-1:length(R1_marker_list_f3.5)
print("Constructing Linking Graph!")
BCV_new_stat<-BCV_1fold_test_fast(marker_list=R1_marker_list_f3.5,data_c=data_c,pp=BCVpp)
BCV_new_stat_diag<-BCV_1fold_test_fast_diag(marker_list=R1_marker_list_f3.5,data_c=data_c)

tg_d<-which(apply(BCV_new_stat_diag[[1]]+BCV_new_stat_diag[[2]]+BCV_new_stat_diag[[3]],1,sum)>3)
BCV_new_stat_diff<-BCV_new_stat
for(kk in 1:length(BCV_new_stat_diff))
{
	BCV_new_stat_diff[[kk]]<-BCV_new_stat_diff[[kk]]*0
	for(i in 1:nrow(BCV_new_stat[[kk]]))
	{
		for(j in 1:ncol(BCV_new_stat[[kk]]))
		{
			if(i<j)
			{
				BCV_new_stat_diff[[kk]][j,i]<-BCV_new_stat_diff[[kk]][i,j]<-BCV_new_stat[[kk]][i,j]-BCV_new_stat_diag[[kk]][i,i]-BCV_new_stat_diag[[kk]][j,j]
			}
		}
	}
}

tg_ids<-which(apply(BCV_new_stat_diff[[1]]<0,1,sum)>(cut_pp*(nrow(BCV_new_stat_diff[[1]]))))
ccc<-BCV_new_stat_diff[[1]]
if(length(tg_ids)>0)
{
	ccc<-BCV_new_stat_diff[[1]][-tg_ids,-tg_ids]
}
ddd<-apply(ccc>1,1,sum)
ddd0<-ddd[order(-ddd)]
eee<-ccc
for(i in 1:length(ddd0))
{
	cc<-names(ddd0)[i]
	if(apply(eee>0,1,sum)[cc]>(cut_pp*nrow(eee)))
	{
		eee<-eee[setdiff(rownames(eee),cc),setdiff(colnames(eee),cc)]	
	}
}

colors = c(-100:100)/50
my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n =200)

ccc<-sign((((BCV_new_stat[[1]][rownames(eee),colnames(eee)]==1)+
(BCV_new_stat[[3]][rownames(eee),colnames(eee)]==1)+
(BCV_new_stat[[2]][rownames(eee),colnames(eee)]==1))==3)+ccc0[rownames(eee),colnames(eee)])
diag(ccc)<-1
heatmap.2(ccc
,Rowv=T,Colv =T,scale="none",main="",
        col=my_palette,breaks=colors,density.info="none",dendrogram="both",
        trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

R1_marker_list_f3.5_0<-R1_marker_list_f3.5
tg_ids<-as.numeric(as.character(rownames(ccc)))
R1_marker_list_f3.5_0<-list()
for(i in 1:length(tg_ids))
{
	R1_marker_list_f3.5_0[[i]]<-R1_marker_list_f3.5[[tg_ids[i]]]
}
marker_list<-R1_marker_list_f3.5_0
graph_links<-ccc
	g1 <- graph_from_adjacency_matrix(ccc)
	largest_cliques(g1)
	cg<-clusters(g1)
	tg_ccc<-as.numeric(names(which(cg$membership==5)))
	tg_clusters_all<-unique(cg$membership)
	cluster_number<-c()
	for(i in 1:length(tg_clusters_all))
	{
       	 cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
	}
	R1_marker_merged_gc<-list()
	for(i in 1:length(tg_clusters_all))
	{
        	ccc<-c()
        	tgs<-which(cg$membership==tg_clusters_all[i])
        	for(j in 1:length(tgs))
       	{
                ccc<-c(ccc,marker_list[[tgs[j]]])
        	}
        	R1_marker_merged_gc[[i]]<-unique(ccc)
	}
	tg1c<-which(cluster_number==1)
	tgoc<-which(cluster_number>1)
	oc_1c_cor<-c()
	if((length(tg1c)>0)&(length(tgoc)>0))
	  {
        oc_1c_cor<-matrix(0,length(tg1c),length(tgoc))
        rownames(oc_1c_cor)<-tg1c
        colnames(oc_1c_cor)<-tgoc
        for(i in 1:length(tg1c))
        {
        for(j in 1:length(tgoc))
        {
                tg1<-R1_marker_merged_gc[[tg1c[i]]]
                tg2<-R1_marker_merged_gc[[tgoc[j]]]
                k1<-sum(BCV_ttest(data_c[tg2,],rounds=40)<0.001)
                k2<-sum(BCV_ttest(data_c[tg1,],rounds=40)<0.001)
                k3<-sum(BCV_ttest(data_c[union(tg1,tg2),],rounds=40)<0.001)
                cc1<-apply(immune_cell_uni_table[tg2,],2,mean)
                cc2<-apply(immune_cell_uni_table[tg1,],2,mean)
                dd1<-c(cc2[which(cc2==max(cc2))[1]]-cc1[which(cc2==max(cc2))[1]])
                dd2<-c(cc1[which(cc1==max(cc1))[1]]-cc2[which(cc1==max(cc1))[1]])
                #print(c(i,j))
                #print(c(dd1,dd2))
                if((k3<=k1)&(dd1<0.5)&(dd2<0.5))
                {
                        oc_1c_cor[i,j]<-1
                }
        }
        }
	  }
Linking_graph_info<-list(graph_links,marker_list,oc_1c_cor,R1_marker_list_f2.5,R1_marker_list_f3.5)
names(Linking_graph_info)<-c("graph_links","marker_list","oc_1c_cor","R1_marker_list_f2.5","R1_marker_list_f3.5")
return(Linking_graph_info)
}


R1_list_filtering_step25_build_LG2<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3,BCVpp=1e-5)
{
data_CORS_cancer<-data_c
module_info<-R1_filter_step1_results[[4]]
module_rank<-R1_filter_step1_results[[5]]
R1_marker_list_f2.5<-clean_rank1_module(data_c,module_info,module_rank,st0=4)
cut1<-cut10
ccc<-compute_min_jaccard(R1_marker_list_f2.5)
ccc0<-ccc>cut1
names(R1_marker_list_f2.5)<-1:length(R1_marker_list_f2.5)
print("Run 1st round BCV and Jaccard distance for merging bases!")
module_merge_info<-compute_BCV_1fold_test(R1_marker_list_f2.5,data_c,ccc0)
tg_cluster_genes<-module_merge_info[[2]]
pp_all<-c()
for(i in 1:length(tg_cluster_genes))
{
        pp<-sum(BCV_ttest(data_c[tg_cluster_genes[[i]],],maxrank0=20)<0.001)  
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-pp_all
R1_marker_list_f3.5<-clean_rank1_module(data_c,tg_cluster_genes,pp_R1_marker_list_f3,st0=4)
cut1=cut10
ccc<-compute_min_jaccard(R1_marker_list_f3.5)
ccc0<-ccc>cut1
rownames(ccc0)<-1:nrow(ccc0)
colnames(ccc0)<-1:ncol(ccc0)

names(R1_marker_list_f3.5)<-1:length(R1_marker_list_f3.5)
print("Constructing Linking Graph!")
BCV_new_stat<-BCV_1fold_test_fast(marker_list=R1_marker_list_f3.5,data_c=data_c,pp=BCVpp)
BCV_new_stat_diag<-BCV_1fold_test_fast_diag(marker_list=R1_marker_list_f3.5,data_c=data_c)

tg_d<-which(apply(BCV_new_stat_diag[[1]]+BCV_new_stat_diag[[2]]+BCV_new_stat_diag[[3]],1,sum)>3)
BCV_new_stat_diff<-BCV_new_stat
for(kk in 1:length(BCV_new_stat_diff))
{
	BCV_new_stat_diff[[kk]]<-BCV_new_stat_diff[[kk]]*0
	for(i in 1:nrow(BCV_new_stat[[kk]]))
	{
		for(j in 1:ncol(BCV_new_stat[[kk]]))
		{
			if(i<j)
			{
				BCV_new_stat_diff[[kk]][j,i]<-BCV_new_stat_diff[[kk]][i,j]<-BCV_new_stat[[kk]][i,j]-BCV_new_stat_diag[[kk]][i,i]-BCV_new_stat_diag[[kk]][j,j]
			}
		}
	}
}

tg_ids<-which(apply(BCV_new_stat_diff[[1]]<0,1,sum)>(cut_pp*(nrow(BCV_new_stat_diff[[1]]))))
ccc<-BCV_new_stat_diff[[1]]
if(length(tg_ids)>0)
{
	ccc<-BCV_new_stat_diff[[1]][-tg_ids,-tg_ids]
}
ddd<-apply(ccc>1,1,sum)
ddd0<-ddd[order(-ddd)]
eee<-ccc
for(i in 1:length(ddd0))
{
	cc<-names(ddd0)[i]
	if(apply(eee>0,1,sum)[cc]>(cut_pp*nrow(eee)))
	{
		eee<-eee[setdiff(rownames(eee),cc),setdiff(colnames(eee),cc)]	
	}
}

colors = c(-100:100)/50
my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n =200)

ccc<-sign((BCV_new_stat[[1]][rownames(eee),colnames(eee)]==1)+
ccc0[rownames(eee),colnames(eee)])
diag(ccc)<-1
heatmap.2(ccc
,Rowv=T,Colv =T,scale="none",main="",
        col=my_palette,breaks=colors,density.info="none",dendrogram="both",
        trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)

R1_marker_list_f3.5_0<-R1_marker_list_f3.5
tg_ids<-as.numeric(as.character(rownames(ccc)))
R1_marker_list_f3.5_0<-list()
for(i in 1:length(tg_ids))
{
	R1_marker_list_f3.5_0[[i]]<-R1_marker_list_f3.5[[tg_ids[i]]]
}
marker_list<-R1_marker_list_f3.5_0
graph_links<-ccc
	g1 <- graph_from_adjacency_matrix(ccc)
	largest_cliques(g1)
	cg<-clusters(g1)
	tg_ccc<-as.numeric(names(which(cg$membership==5)))
	tg_clusters_all<-unique(cg$membership)
	cluster_number<-c()
	for(i in 1:length(tg_clusters_all))
	{
       	 cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
	}
	R1_marker_merged_gc<-list()
	for(i in 1:length(tg_clusters_all))
	{
        	ccc<-c()
        	tgs<-which(cg$membership==tg_clusters_all[i])
        	for(j in 1:length(tgs))
       	{
                ccc<-c(ccc,marker_list[[tgs[j]]])
        	}
        	R1_marker_merged_gc[[i]]<-unique(ccc)
	}
	tg1c<-which(cluster_number==1)
	tgoc<-which(cluster_number>1)
	oc_1c_cor<-c()
	  if((length(tg1c)>0)&(length(tgoc)>0))
	  {
        oc_1c_cor<-matrix(0,length(tg1c),length(tgoc))
        rownames(oc_1c_cor)<-tg1c
        colnames(oc_1c_cor)<-tgoc
        for(i in 1:length(tg1c))
        {
        for(j in 1:length(tgoc))
        {
                tg1<-R1_marker_merged_gc[[tg1c[i]]]
                tg2<-R1_marker_merged_gc[[tgoc[j]]]
                k1<-sum(BCV_ttest(data_c[tg2,],rounds=40)<0.001)
                k2<-sum(BCV_ttest(data_c[tg1,],rounds=40)<0.001)
                k3<-sum(BCV_ttest(data_c[union(tg1,tg2),],rounds=40)<0.001)
                cc1<-apply(immune_cell_uni_table[tg2,],2,mean)
                cc2<-apply(immune_cell_uni_table[tg1,],2,mean)
                dd1<-c(cc2[which(cc2==max(cc2))[1]]-cc1[which(cc2==max(cc2))[1]])
                dd2<-c(cc1[which(cc1==max(cc1))[1]]-cc2[which(cc1==max(cc1))[1]])
                #print(c(i,j))
                #print(c(dd1,dd2))
                if((k3<=k1)&(dd1<0.5)&(dd2<0.5))
                {
                        oc_1c_cor[i,j]<-1
                }
        }
        }
	  }

Linking_graph_info<-list(graph_links,marker_list,oc_1c_cor,R1_marker_list_f2.5,R1_marker_list_f3.5)
names(Linking_graph_info)<-c("graph_links","marker_list","oc_1c_cor","R1_marker_list_f2.5","R1_marker_list_f3.5")
return(Linking_graph_info)
}


R1_list_filtering_step35_buildNMF_pure1<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
{
	graph_links<-ccc<-Linking_graph_info[[1]]
	marker_list<-Linking_graph_info[[2]]
	g1 <- graph_from_adjacency_matrix(ccc)
	largest_cliques(g1)
	cg<-clusters(g1)
	tg_ccc<-as.numeric(names(which(cg$membership==5)))
	tg_clusters_all<-unique(cg$membership)
	cluster_number<-c()
	for(i in 1:length(tg_clusters_all))
	{
       	 cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
	}
	R1_marker_merged_gc<-list()
	for(i in 1:length(tg_clusters_all))
	{
        	ccc<-c()
        	tgs<-which(cg$membership==tg_clusters_all[i])
        	for(j in 1:length(tgs))
       	{
                ccc<-c(ccc,marker_list[[tgs[j]]])
        	}
        	R1_marker_merged_gc[[i]]<-unique(ccc)
	}
	tg1c<-which(cluster_number==1)
	tgoc<-which(cluster_number>1)
	oc_1c_cor<-Linking_graph_info[[3]]
	tg1c_new<-c()
	tg1c_merge<-c()
	if(length(oc_1c_cor)>0)
	{
		if(sum(apply(oc_1c_cor,1,sum)==0)>0)
		{
        		tg1c_new<-tg1c[which(apply(oc_1c_cor,1,sum)==0)]
		}
		tg1c_merge<-setdiff(tg1c,tg1c_new)
	}
	if(length(oc_1c_cor)==0)
	{
		tg1c_new<-tg1c
	}
	All_genes<-c()
	for(i in 1:length(marker_list))
	{
        	All_genes<-c(All_genes,marker_list[[i]])
	}
	All_genes<-sort(unique(All_genes))

NMF_indi_all<-c()
LG_id_all<-list()
NMF_indi_list<-list()
N<-0
print("Working on Linking graph for NMF constraints!")
if(length(tgoc)>0)
{
for(ii in 1:length(tgoc))
{
	N<-N+1
        print(c(ii,length(tgoc)))
        tg1c_m_c<-tg1c[which(oc_1c_cor[,ii]==1)]
        tgs<-which(cg$membership==tgoc[ii])
	  tgs0<-tgs
        if(length(tg1c_m_c)>0)
        {
                for(j in 1:length(tg1c_m_c))
                {
                        tgs<-c(tgs,which(cg$membership==tg1c_m_c[j]))
                }
        }
        LG_id_all[[N]]<-list(tgs0,tgs)
	  tg_M_all0<-c()
        for(j in 1:length(tgs0))
        {
                tg_M_all0<-c(tg_M_all0,marker_list[[tgs0[j]]])
        }
        tg_M_all0<-unique(tg_M_all0)
        tg_M_all<-c()
        for(j in 1:length(tgs))
        {
                tg_M_all<-c(tg_M_all,marker_list[[tgs[j]]])
        }
        tg_M_all<-unique(tg_M_all)
        pp<-c()
        for(k in 1:21)
        {
                pp<-c(pp,sum(BCV_ttest(data_c[tg_M_all,],rounds=40,maxrank0=20)<0.01))
        }
        pp0<-floor(median(pp))
        Base_all<-c()
        for(j in 1:length(tgs))
        {
                Base_all<-rbind(Base_all,apply(data_c[marker_list[[tgs[j]]],],2,mean))
        }
        NMF_label_c<-c()
        if(pp0<=1)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste(tgoc[ii],"All1_0",N,sep="_")
        }
        if(pp0>1)
        {
        stat_all<-rep(1,length(tgs))
        names(stat_all)<-tgs
        if(length(tgs)>=3)
        {
        ddd<-combn(length(tgs),3)
        comp_ids<-c()
        stat_all<-rep(0,length(tgs))
        names(stat_all)<-tgs
        for(i in 1:ncol(ddd))
        {
                tg_genes<-unique(c(marker_list[[tgs[ddd[1,i]]]],marker_list[[tgs[ddd[2,i]]]],marker_list[[tgs[ddd[3,i]]]]))
                pp<-sum(BCV_ttest(data_c[tg_genes,],rounds=20,maxrank0=5)<0.01)
                if(pp==2)
                {
                        comp_ids<-c(comp_ids,i)
                        #print(i)
                        lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
                        if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[2]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==-1))
                        {
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                }
        }
        }
	  stat_all<-stat_all#[as.character(tgs0)]
        ccc<-stat_all[order(-stat_all)]
        ccc<-ccc[which(ccc!=0)]
        specific_cell_marker_ids<-c()
        if(length(ccc)>0)
        {
        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[1]))
        for(j in 1:length(ccc))
        {
                if(sum(graph_links[specific_cell_marker_ids,as.numeric(names(ccc)[j])])==0)
                {
                        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[j]))
                }
        }
        if(length(specific_cell_marker_ids)>pp0)
        {
                specific_cell_marker_ids<-specific_cell_marker_ids[1:pp0]
        }
        NMF_label_c<-c()
        if(length(specific_cell_marker_ids)>1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all
                colnames(NMF_label_c)<-paste("Selected",specific_cell_marker_ids,sep="_")
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[marker_list[[specific_cell_marker_ids[i]]],i]<-1
                }
                tg_all_genes<-setdiff(tg_M_all,unique(tg_specific_gene))
                NMF_label_c[tg_all_genes,]<-NMF_label_c[tg_all_genes,]+0
        }
        if(length(specific_cell_marker_ids)==1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all0),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all0
                colnames(NMF_label_c)<-specific_cell_marker_ids
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[intersect(rownames(NMF_label_c),marker_list[[specific_cell_marker_ids[i]]]),i]<-1
                }
                NMF_label_c<-cbind(NMF_label_c,rep(1,length(tg_M_all0)))
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_1",N,sep="_")
        }
        }
        if((length(specific_cell_marker_ids)==0)&pp<3)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_2",N,sep="_")
        }
        if((length(specific_cell_marker_ids)==0)&pp>=3)
        {
                NMF_label_c<-c()
                NMF_label_c0<-rep(0,length(tg_M_all))
                names(NMF_label_c0)<-tg_M_all
                for(j in 1:length(tgs))
                {
                        NMF_label_ccc<-NMF_label_c0
                        NMF_label_ccc[marker_list[[tgs[j]]]]<-1
                        NMF_label_c<-cbind(NMF_label_c,NMF_label_ccc)
                }
                colnames(NMF_label_c)<-paste("tgoc",tgs,"ss1",N,sep="_")
        }
        } 
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
	  if(ncol(NMF_label_c)>1)
	  {
		NMF_label_c<-NMF_label_c[which(apply(NMF_label_c,1,sum)>0),]
        }
	  NMF_indi_list[[N]]<-NMF_label_c
}
}
if(length(tg1c_new)>0)
{
for(ii in 1:length(tg1c_new))
{
	  N<-N+1
        tgs<-which(cg$membership==tg1c_new[ii])
	  LG_id_all[[N]]<-list(tgs,tgs)
        tg_M_all<-marker_list[[tgs]]
        NMF_label_c<-rep(1,length(tg_M_all))
        names(NMF_label_c)<-tg_M_all
        NMF_label_c<-as.matrix(NMF_label_c)
        colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tg1c",tg1c_new[ii],"All1",sep="_")
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
	  NMF_indi_list[[N]]<-NMF_label_c
}
}
NMF_indi_all_f<-c()
ccc<-unique(colnames(NMF_indi_all))
for(i in 1:length(ccc))
{
	tg_ids<-which(colnames(NMF_indi_all)==ccc[i])
	if(length(tg_ids)==1)
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,NMF_indi_all[,tg_ids])
	}
	else
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,sign(apply(NMF_indi_all[,tg_ids],1,sum))*1)
	}
}
colnames(NMF_indi_all_f)<-ccc
to_d<-c()
for(i in 1:ncol(NMF_indi_all_f))
{
	for(j in 1:ncol(NMF_indi_all_f))
	{
		if(i<j)
		{
			cc<-sum(abs(NMF_indi_all_f[,i]-NMF_indi_all_f[,j]))
			if(cc==0)
			{
				to_d<-c(to_d,j)
			}
		}
	}
}
if(length(to_d)>0)
{
	NMF_indi_all_f<-NMF_indi_all_f[,-to_d]
}
NMF_indi_all_f<-NMF_indi_all_f[which(apply(NMF_indi_all_f,1,sum)>0),]
R1_filter_step3_results<-list(NMF_indi_all_f,NMF_indi_list,LG_id_all,graph_links,g1,cg)
names(R1_filter_step3_results)<-c("NMF_indi_all","NMF_indi_list","LG_id_all","graph_links","constructed graph","graph cluster")
return(R1_filter_step3_results)
}

R1_list_filtering_step35_buildNMF_pure2<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
{
	graph_links<-ccc<-Linking_graph_info[[1]]
	marker_list<-Linking_graph_info[[2]]
	g1 <- graph_from_adjacency_matrix(ccc)
	largest_cliques(g1)
	cg<-clusters(g1)
	tg_ccc<-as.numeric(names(which(cg$membership==5)))
	tg_clusters_all<-unique(cg$membership)
	cluster_number<-c()
	for(i in 1:length(tg_clusters_all))
	{
       	 cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
	}
	R1_marker_merged_gc<-list()
	for(i in 1:length(tg_clusters_all))
	{
        	ccc<-c()
        	tgs<-which(cg$membership==tg_clusters_all[i])
        	for(j in 1:length(tgs))
       	{
                ccc<-c(ccc,marker_list[[tgs[j]]])
        	}
        	R1_marker_merged_gc[[i]]<-unique(ccc)
	}
	tg1c<-which(cluster_number==1)
	tgoc<-which(cluster_number>1)
	oc_1c_cor<-Linking_graph_info[[3]]
	tg1c_new<-c()
	tg1c_merge<-c()
	if(length(oc_1c_cor)>0)
	{
	if(sum(apply(oc_1c_cor,1,sum)==0)>0)
	{
        	tg1c_new<-tg1c[which(apply(oc_1c_cor,1,sum)==0)]
	}
	tg1c_merge<-setdiff(tg1c,tg1c_new)
	}
	if(length(oc_1c_cor)==0)
	{
		tg1c_new<-tg1c
	}
	All_genes<-c()
	for(i in 1:length(marker_list))
	{
        	All_genes<-c(All_genes,marker_list[[i]])
	}
	All_genes<-sort(unique(All_genes))

NMF_indi_all<-c()
LG_id_all<-list()
NMF_indi_list<-list()
N<-0
print("Working on Linking graph for NMF constraints!")
if(length(tgoc)>0)
{
for(ii in 1:length(tgoc))
{
	N<-N+1
        print(c(ii,length(tgoc)))
        tg1c_m_c<-tg1c[which(oc_1c_cor[,ii]==1)]
        tgs<-which(cg$membership==tgoc[ii])
	  tgs0<-tgs
        if(length(tg1c_m_c)>0)
        {
                for(j in 1:length(tg1c_m_c))
                {
                        tgs<-c(tgs,which(cg$membership==tg1c_m_c[j]))
                }
        }
        LG_id_all[[N]]<-list(tgs0,tgs)
	  tg_M_all0<-c()
        for(j in 1:length(tgs0))
        {
                tg_M_all0<-c(tg_M_all0,marker_list[[tgs0[j]]])
        }
        tg_M_all0<-unique(tg_M_all0)
        tg_M_all<-c()
        for(j in 1:length(tgs))
        {
                tg_M_all<-c(tg_M_all,marker_list[[tgs[j]]])
        }
        tg_M_all<-unique(tg_M_all)
        pp<-c()
        for(k in 1:21)
        {
                pp<-c(pp,sum(BCV_ttest(data_c[tg_M_all,],rounds=40,maxrank0=20)<0.01))
        }
        pp0<-floor(median(pp))
        Base_all<-c()
        for(j in 1:length(tgs))
        {
                Base_all<-rbind(Base_all,apply(data_c[marker_list[[tgs[j]]],],2,mean))
        }
        NMF_label_c<-c()
        if(pp0<=1)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste(tgoc[ii],"All1_0",N,sep="_")
        }
        if(pp0>1)
        {
        stat_all<-rep(1,length(tgs))
        names(stat_all)<-tgs
        if(length(tgs)>=3)
        {
        ddd<-combn(length(tgs),3)
        comp_ids<-c()
        stat_all<-rep(0,length(tgs))
        names(stat_all)<-tgs
        for(i in 1:ncol(ddd))
        {
                tg_genes<-unique(c(marker_list[[tgs[ddd[1,i]]]],marker_list[[tgs[ddd[2,i]]]],marker_list[[tgs[ddd[3,i]]]]))
                pp<-sum(BCV_ttest(data_c[tg_genes,],rounds=20,maxrank0=5)<0.01)
                if(pp==2)
                {
                        comp_ids<-c(comp_ids,i)
                        #print(i)
                        lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
                        if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                }
        }
        }
	  stat_all<-stat_all#[as.character(tgs0)]
        ccc<-stat_all[order(-stat_all)]
        ccc<-ccc[which(ccc!=0)]
        specific_cell_marker_ids<-c()
        if(length(ccc)>0)
        {
        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[1]))
        for(j in 1:length(ccc))
        {
                if(sum(graph_links[specific_cell_marker_ids,as.numeric(names(ccc)[j])])==0)
                {
                        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[j]))
                }
        }
        if(length(specific_cell_marker_ids)>pp0)
        {
                specific_cell_marker_ids<-specific_cell_marker_ids[1:pp0]
        }
        NMF_label_c<-c()
        if(length(specific_cell_marker_ids)>1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all
                colnames(NMF_label_c)<-paste("Selected",specific_cell_marker_ids,sep="_")
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[marker_list[[specific_cell_marker_ids[i]]],i]<-1
                }
                tg_all_genes<-setdiff(tg_M_all,unique(tg_specific_gene))
                NMF_label_c[tg_all_genes,]<-NMF_label_c[tg_all_genes,]+1
        }
        if(length(specific_cell_marker_ids)==1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all0),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all0
                colnames(NMF_label_c)<-specific_cell_marker_ids
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[intersect(rownames(NMF_label_c),marker_list[[specific_cell_marker_ids[i]]]),i]<-1
                }
                NMF_label_c<-cbind(NMF_label_c,rep(1,length(tg_M_all0)))
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_1",N,sep="_")
        }
        }
        if((length(specific_cell_marker_ids)==0)&pp<3)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_2",N,sep="_")
        }
        if((length(specific_cell_marker_ids)==0)&pp>=3)
        {
                NMF_label_c<-c()
                NMF_label_c0<-rep(0,length(tg_M_all))
                names(NMF_label_c0)<-tg_M_all
                for(j in 1:length(tgs))
                {
                        NMF_label_ccc<-NMF_label_c0
                        NMF_label_ccc[marker_list[[tgs[j]]]]<-1
                        NMF_label_c<-cbind(NMF_label_c,NMF_label_ccc)
                }
                colnames(NMF_label_c)<-paste("tgoc",tgs,"ss1",N,sep="_")
        }
        } 
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
	  if(ncol(NMF_label_c)>1)
	  {
		NMF_label_c<-NMF_label_c[which(apply(NMF_label_c,1,sum)>0),]
        }
	  NMF_indi_list[[N]]<-NMF_label_c
}
}
if(length(tg1c_new)>0)
{
for(ii in 1:length(tg1c_new))
{
	  N<-N+1
        tgs<-which(cg$membership==tg1c_new[ii])
	  LG_id_all[[N]]<-list(tgs,tgs)
        tg_M_all<-marker_list[[tgs]]
        NMF_label_c<-rep(1,length(tg_M_all))
        names(NMF_label_c)<-tg_M_all
        NMF_label_c<-as.matrix(NMF_label_c)
        colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tg1c",tg1c_new[ii],"All1",sep="_")
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
	  NMF_indi_list[[N]]<-NMF_label_c
}
}
NMF_indi_all_f<-c()
ccc<-unique(colnames(NMF_indi_all))
for(i in 1:length(ccc))
{
	tg_ids<-which(colnames(NMF_indi_all)==ccc[i])
	if(length(tg_ids)==1)
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,NMF_indi_all[,tg_ids])
	}
	else
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,sign(apply(NMF_indi_all[,tg_ids],1,sum))*1)
	}
}
colnames(NMF_indi_all_f)<-ccc
to_d<-c()
for(i in 1:ncol(NMF_indi_all_f))
{
	for(j in 1:ncol(NMF_indi_all_f))
	{
		if(i<j)
		{
			cc<-sum(abs(NMF_indi_all_f[,i]-NMF_indi_all_f[,j]))
			if(cc==0)
			{
				to_d<-c(to_d,j)
			}
		}
	}
}
if(length(to_d)>0)
{
	NMF_indi_all_f<-NMF_indi_all_f[,-to_d]
}
NMF_indi_all_f<-NMF_indi_all_f[which(apply(NMF_indi_all_f,1,sum)>0),]
R1_filter_step3_results<-list(NMF_indi_all_f,NMF_indi_list,LG_id_all,graph_links,g1,cg)
names(R1_filter_step3_results)<-c("NMF_indi_all","NMF_indi_list","LG_id_all","graph_links","constructed graph","graph cluster")
return(R1_filter_step3_results)
}


generate_oc_1c_cor<-function(Linking_graph_info,data_c=data_CORS_cancer,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
{
	graph_links<-ccc<-Linking_graph_info[[1]]
	marker_list<-Linking_graph_info[[2]]
	g1 <- graph_from_adjacency_matrix(ccc)
	largest_cliques(g1)
	cg<-clusters(g1)
	tg_ccc<-as.numeric(names(which(cg$membership==5)))
	tg_clusters_all<-unique(cg$membership)
	cluster_number<-c()
	for(i in 1:length(tg_clusters_all))
	{
       	 cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
	}
	R1_marker_merged_gc<-list()
	for(i in 1:length(tg_clusters_all))
	{
        	ccc<-c()
        	tgs<-which(cg$membership==tg_clusters_all[i])
        	for(j in 1:length(tgs))
       	{
                ccc<-c(ccc,marker_list[[tgs[j]]])
        	}
        	R1_marker_merged_gc[[i]]<-unique(ccc)
	}
	tg1c<-which(cluster_number==1)
	tgoc<-which(cluster_number>1)
	oc_1c_cor<-matrix(0,length(tg1c),length(tgoc))
	rownames(oc_1c_cor)<-tg1c
	colnames(oc_1c_cor)<-tgoc
	for(i in 1:length(tg1c))
	{
        for(j in 1:length(tgoc))
        {
                tg1<-R1_marker_merged_gc[[tg1c[i]]]
                tg2<-R1_marker_merged_gc[[tgoc[j]]]
                k1<-sum(BCV_ttest(data_c[tg2,],rounds=40)<0.001)
                k2<-sum(BCV_ttest(data_c[tg1,],rounds=40)<0.001)
                k3<-sum(BCV_ttest(data_c[union(tg1,tg2),],rounds=40)<0.001)
                cc1<-apply(immune_cell_uni_table[tg2,],2,mean)
                cc2<-apply(immune_cell_uni_table[tg1,],2,mean)
                dd1<-c(cc2[which(cc2==max(cc2))[1]]-cc1[which(cc2==max(cc2))[1]])
                dd2<-c(cc1[which(cc1==max(cc1))[1]]-cc2[which(cc1==max(cc1))[1]])
                #print(c(i,j))
                #print(c(dd1,dd2))
                if((k3<=k1)&(dd1<0.5)&(dd2<0.5))
                {
                        oc_1c_cor[i,j]<-1
                }
        }
	}
	return(oc_1c_cor)
}

R1_list_filtering_step3_buildNMF_new<-function(R1_filter_step1_results,data_c=data_CORS_cancer,cut10=0.8,immune_cell_uni_table=immune_cell_uni_table0_GPL570,cut_pp=0.3)
{
data_CORS_cancer<-data_c
module_info<-R1_filter_step1_results[[4]]
module_rank<-R1_filter_step1_results[[5]]
R1_marker_list_f2.5<-clean_rank1_module(data_c,module_info,module_rank,st0=4)
cut1<-cut10
ccc<-compute_min_jaccard(R1_marker_list_f2.5)
ccc0<-ccc>cut1
names(R1_marker_list_f2.5)<-1:length(R1_marker_list_f2.5)
print("Run 1st round BCV and Jaccard distance for merging bases!")
module_merge_info<-compute_BCV_1fold_test(R1_marker_list_f2.5,data_c,ccc0)
tg_cluster_genes<-module_merge_info[[2]]
pp_all<-c()
for(i in 1:length(tg_cluster_genes))
{
        pp<-sum(BCV_ttest(data_c[tg_cluster_genes[[i]],],maxrank0=20)<0.001)  
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-pp_all
R1_marker_list_f3.5<-clean_rank1_module(data_c,tg_cluster_genes,pp_R1_marker_list_f3,st0=4)
cut1=cut10
ccc<-compute_min_jaccard(R1_marker_list_f3.5)
ccc0<-ccc>cut1
rownames(ccc0)<-1:nrow(ccc0)
colnames(ccc0)<-1:ncol(ccc0)

names(R1_marker_list_f3.5)<-1:length(R1_marker_list_f3.5)
print("Constructing Linking Graph!")
BCV_new_stat<-BCV_1fold_test_fast(marker_list=R1_marker_list_f3.5,data_c=data_c)
BCV_new_stat_diag<-BCV_1fold_test_fast_diag(marker_list=R1_marker_list_f3.5,data_c=data_c)

tg_d<-which(apply(BCV_new_stat_diag[[1]]+BCV_new_stat_diag[[2]]+BCV_new_stat_diag[[3]],1,sum)>3)
BCV_new_stat_diff<-BCV_new_stat
for(kk in 1:length(BCV_new_stat_diff))
{
	BCV_new_stat_diff[[kk]]<-BCV_new_stat_diff[[kk]]*0
	for(i in 1:nrow(BCV_new_stat[[kk]]))
	{
		for(j in 1:ncol(BCV_new_stat[[kk]]))
		{
			if(i<j)
			{
				BCV_new_stat_diff[[kk]][j,i]<-BCV_new_stat_diff[[kk]][i,j]<-BCV_new_stat[[kk]][i,j]-BCV_new_stat_diag[[kk]][i,i]-BCV_new_stat_diag[[kk]][j,j]
			}
		}
	}
}
tg_ids<-which(apply(BCV_new_stat_diff[[1]]<0,1,sum)>(cut_pp*(nrow(BCV_new_stat_diff[[1]]))))
ccc<-BCV_new_stat_diff[[1]]
if(length(tg_ids)>0)
{
	ccc<-BCV_new_stat_diff[[1]][-tg_ids,-tg_ids]
}
ddd<-apply(ccc>1,1,sum)
ddd0<-ddd[order(-ddd)]
eee<-ccc
for(i in 1:length(ddd0))
{
	cc<-names(ddd0)[i]
	if(apply(eee>0,1,sum)[cc]>(cut_pp*nrow(eee)))
	{
		eee<-eee[setdiff(rownames(eee),cc),setdiff(colnames(eee),cc)]	
	}
}

colors = c(-100:100)/50
my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n =200)

ccc<-sign((BCV_new_stat[[1]][rownames(eee),colnames(eee)]==1)+
(BCV_new_stat[[3]][rownames(eee),colnames(eee)]==1)+
(BCV_new_stat[[2]][rownames(eee),colnames(eee)]==1)+ccc0[rownames(eee),colnames(eee)])

diag(ccc)<-1
heatmap.2(ccc
,Rowv=T,Colv =T,scale="none",main="",
        col=my_palette,breaks=colors,density.info="none",dendrogram="both",
        trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)


R1_marker_list_f3.5_0<-R1_marker_list_f3.5
tg_ids<-as.numeric(as.character(rownames(ccc)))
R1_marker_list_f3.5_0<-list()
for(i in 1:length(tg_ids))
{
	R1_marker_list_f3.5_0[[i]]<-R1_marker_list_f3.5[[tg_ids[i]]]
}
marker_list<-R1_marker_list_f3.5_0
graph_links<-ccc

g1 <- graph_from_adjacency_matrix(ccc)
largest_cliques(g1)
cg<-clusters(g1)
tg_ccc<-as.numeric(names(which(cg$membership==5)))
tg_clusters_all<-unique(cg$membership)
cluster_number<-c()
for(i in 1:length(tg_clusters_all))
{
        cluster_number<-c(cluster_number,sum(cg$membership==tg_clusters_all[i]))
}
R1_marker_merged_gc<-list()
for(i in 1:length(tg_clusters_all))
{
        ccc<-c()
        tgs<-which(cg$membership==tg_clusters_all[i])
        for(j in 1:length(tgs))
        {
                ccc<-c(ccc,marker_list[[tgs[j]]])
        }
        R1_marker_merged_gc[[i]]<-unique(ccc)
}
tg1c<-which(cluster_number==1)
tgoc<-which(cluster_number>1)
oc_1c_cor<-matrix(0,length(tg1c),length(tgoc))
rownames(oc_1c_cor)<-tg1c
colnames(oc_1c_cor)<-tgoc
for(i in 1:length(tg1c))
{
        for(j in 1:length(tgoc))
        {
                tg1<-R1_marker_merged_gc[[tg1c[i]]]
                tg2<-R1_marker_merged_gc[[tgoc[j]]]
                k1<-sum(BCV_ttest(data_c[tg2,],rounds=40)<0.001)
                k2<-sum(BCV_ttest(data_c[tg1,],rounds=40)<0.001)
                k3<-sum(BCV_ttest(data_c[union(tg1,tg2),],rounds=40)<0.001)
                cc1<-apply(immune_cell_uni_table[tg2,],2,mean)
                cc2<-apply(immune_cell_uni_table[tg1,],2,mean)
                dd1<-c(cc2[which(cc2==max(cc2))[1]]-cc1[which(cc2==max(cc2))[1]])
                dd2<-c(cc1[which(cc1==max(cc1))[1]]-cc2[which(cc1==max(cc1))[1]])
                #print(c(i,j))
                #print(c(dd1,dd2))
                if((k3<=k1)&(dd1<0.5)&(dd2<0.5))
                {
                        oc_1c_cor[i,j]<-1
                }
        }
}
tg1c_new<-c()
if(sum(apply(oc_1c_cor,1,sum)==0)>0)
{
        tg1c_new<-tg1c[which(apply(oc_1c_cor,1,sum)==0)]
}
tg1c_merge<-setdiff(tg1c,tg1c_new)
All_genes<-c()
for(i in 1:length(marker_list))
{
        All_genes<-c(All_genes,marker_list[[i]])
}
All_genes<-sort(unique(All_genes))

NMF_indi_all<-c()
print("Working on Linking graph for NMF constraints!")
for(ii in 1:length(tgoc))
{
        print(c(ii,length(tgoc)))
        tg1c_m_c<-tg1c[which(oc_1c_cor[,ii]==1)]
        tgs<-which(cg$membership==tgoc[ii])
	  tgs0<-tgs
        if(length(tg1c_m_c)>0)
        {
                for(j in 1:length(tg1c_m_c))
                {
                        tgs<-c(tgs,which(cg$membership==tg1c_m_c[j]))
                }
        }
	  tg_M_all0<-c()
        for(j in 1:length(tgs0))
        {
                tg_M_all0<-c(tg_M_all0,marker_list[[tgs0[j]]])
        }
        tg_M_all0<-unique(tg_M_all0)

        tg_M_all<-c()
        for(j in 1:length(tgs))
        {
                tg_M_all<-c(tg_M_all,marker_list[[tgs[j]]])
        }
        tg_M_all<-unique(tg_M_all)
        pp<-c()
        for(k in 1:21)
        {
                pp<-c(pp,sum(BCV_ttest(data_c[tg_M_all,],rounds=40,maxrank0=20)<0.01))
        }
        pp0<-floor(median(pp))
        Base_all<-c()
        for(j in 1:length(tgs))
        {
                Base_all<-rbind(Base_all,apply(data_c[marker_list[[tgs[j]]],],2,mean))
        }
        NMF_label_c<-c()
        if(pp0<=1)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste(tgoc[ii],"All1",sep="_")
        }
        if(pp0>1)
        {
        stat_all<-rep(1,length(tgs))
        names(stat_all)<-tgs
        if(length(tgs)>=3)
        {
        ddd<-combn(length(tgs),3)
        comp_ids<-c()
        stat_all<-rep(0,length(tgs))
        names(stat_all)<-tgs
        for(i in 1:ncol(ddd))
        {
                tg_genes<-unique(c(marker_list[[tgs[ddd[1,i]]]],marker_list[[tgs[ddd[2,i]]]],marker_list[[tgs[ddd[3,i]]]]))
                pp<-sum(BCV_ttest(data_c[tg_genes,],rounds=20,maxrank0=5)<0.01)
                if(pp==2)
                {
                        comp_ids<-c(comp_ids,i)
                        #print(i)
                        lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
                        if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[2,i]]<-stat_all[ddd[2,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                        if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[1]==1))
                        {
                                stat_all[ddd[3,i]]<-stat_all[ddd[3,i]]+1
                                stat_all[ddd[1,i]]<-stat_all[ddd[1,i]]+1
                        }
                }
        }
        }
	  stat_all<-stat_all#[as.character(tgs0)]
        ccc<-stat_all[order(-stat_all)]
        ccc<-ccc[which(ccc!=0)]
        specific_cell_marker_ids<-c()
        if(length(ccc)>0)
        {
        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[1]))
        for(j in 1:length(ccc))
        {
                if(sum(graph_links[specific_cell_marker_ids,as.numeric(names(ccc)[j])])==0)
                {
                        specific_cell_marker_ids<-c(specific_cell_marker_ids,as.numeric(names(ccc)[j]))
                }
        }
        if(length(specific_cell_marker_ids)>pp0)
        {
                specific_cell_marker_ids<-specific_cell_marker_ids[1:pp0]
        }
        NMF_label_c<-c()
        if(length(specific_cell_marker_ids)>1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all
                colnames(NMF_label_c)<-specific_cell_marker_ids
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[marker_list[[specific_cell_marker_ids[i]]],i]<-1
                }
                tg_all_genes<-setdiff(tg_M_all,unique(tg_specific_gene))
                NMF_label_c[tg_all_genes,]<-NMF_label_c[tg_all_genes,]+0
        }
        if(length(specific_cell_marker_ids)==1)
        {
                NMF_label_c<-matrix(0,length(tg_M_all0),length(specific_cell_marker_ids))
                rownames(NMF_label_c)<-tg_M_all0
                colnames(NMF_label_c)<-specific_cell_marker_ids
                tg_specific_gene<-c()
                for(i in 1:length(specific_cell_marker_ids))
                {
                        tg_specific_gene<-c(tg_specific_gene,marker_list[[specific_cell_marker_ids[i]]])
                        NMF_label_c[intersect(rownames(NMF_label_c),marker_list[[specific_cell_marker_ids[i]]]),i]<-1
                }
                NMF_label_c<-cbind(NMF_label_c,rep(1,length(tg_M_all0)))
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_1",sep="_")
        }
        }
        if((length(specific_cell_marker_ids)==0)&pp<3)
        {
                NMF_label_c<-rep(1,length(tg_M_all))
                names(NMF_label_c)<-tg_M_all
                NMF_label_c<-as.matrix(NMF_label_c)
                colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tgoc",tgoc[ii],"All1_2",sep="_")
        }
        if((length(specific_cell_marker_ids)==0)&pp>=3)
        {
                NMF_label_c<-c()
                NMF_label_c0<-rep(0,length(tg_M_all))
                names(NMF_label_c0)<-tg_M_all
                for(j in 1:length(tgs))
                {
                        NMF_label_ccc<-NMF_label_c0
                        NMF_label_ccc[marker_list[[tgs[j]]]]<-1
                        NMF_label_c<-cbind(NMF_label_c,NMF_label_ccc)
                }
                colnames(NMF_label_c)<-paste("tgoc",tgs,"ss1",sep="_")
        }
        }
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
}
if(length(tg1c_new)>0)
{
for(ii in 1:length(tg1c_new))
{
        tgs<-which(cg$membership==tg1c_new[ii])
        tg_M_all<-marker_list[[tgs]]
        NMF_label_c<-rep(1,length(tg_M_all))
        names(NMF_label_c)<-tg_M_all
        NMF_label_c<-as.matrix(NMF_label_c)
        colnames(NMF_label_c)[ncol(NMF_label_c)]<-paste("tg1c",tg1c_new[ii],"All1",sep="_")
        NMF_label_c_all<-matrix(0,length(All_genes),ncol(NMF_label_c))
        colnames(NMF_label_c_all)<-colnames(NMF_label_c)
        rownames(NMF_label_c_all)<-All_genes
        NMF_label_c_all[rownames(NMF_label_c),colnames(NMF_label_c)]<-NMF_label_c
        NMF_indi_all<-cbind(NMF_indi_all,NMF_label_c_all)
}
}
NMF_indi_all_f<-c()
ccc<-unique(colnames(NMF_indi_all))
for(i in 1:length(ccc))
{
	tg_ids<-which(colnames(NMF_indi_all)==ccc[i])
	if(length(tg_ids)==1)
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,NMF_indi_all[,tg_ids])
	}
	else
	{
		NMF_indi_all_f<-cbind(NMF_indi_all_f,sign(apply(NMF_indi_all[,tg_ids],1,sum))*1)
	}
}
colnames(NMF_indi_all_f)<-ccc
NMF_indi_all_f<-NMF_indi_all_f[which(apply(NMF_indi_all_f,1,sum)>0),]
R1_filter_step3_results<-list(NMF_indi_all_f,R1_marker_list_f2.5,R1_marker_list_f3.5,graph_links,g1,cg)
names(R1_filter_step3_results)<-c("NMF_indi_all","R1_marker_list_f2.5","R1_marker_list_f3.5","graph_links","constructed graph","graph cluster")

return(R1_filter_step3_results)
}



Reduce_NMF_redundancy<-function(C_ids,NMF_indi_all,data_CORS_cancer)
{
C_ids_all<-unique(C_ids)
N<-0
Cluster_list<-list()
nn<-c()
for(ii in 1:length(C_ids_all))
{
	C_ids_c<-which(C_ids==C_ids_all[[ii]])
	stat_c<-matrix(0,length(C_ids_c),length(C_ids_c))
	for(i in 1:length(C_ids_c))
	{
		print(i)
		for(j in 1:length(C_ids_c))
		{
			if(i<=j)
			{
			tg_genes1<-names(which(NMF_indi_all[,C_ids_c[i]]==1))
			eee<-cor(t(data_CORS_cancer[tg_genes1,]))
			diag(eee)<-0
			tg_genes1<-names(which(apply(eee,1,mean)>0.2))

			tg_genes2<-names(which(NMF_indi_all[,C_ids_c[j]]==1))
			eee<-cor(t(data_CORS_cancer[tg_genes2,]))
			diag(eee)<-0
			tg_genes2<-names(which(apply(eee,1,mean)>0.2))

			tg_genes<-unique(c(tg_genes1,tg_genes2))
				if(length( tg_genes)>4)
				{
                       	 	stat_c[i,j]<-stat_c[j,i]<-sum(BCV_ttest2(data_CORS_cancer[tg_genes,],rounds=50,msep_cut=0.05)<0.001)
				}
			}
		}
	}
	stat_c0<-(stat_c==1)*1
	rownames(stat_c0)<-C_ids_c
	colnames(stat_c0)<-C_ids_c

	diag(stat_c0)<-1
	graph_links<-stat_c0
	g1 <- graph_from_adjacency_matrix(stat_c0)
	cg<-clusters(g1)
	tg_clusters_all<-unique(cg$membership)
	for(i in 1:length(tg_clusters_all))
	{
		nn<-c(nn,C_ids_all[[ii]])
		N<-N+1
		Cluster_list[[N]]<-as.numeric(names(which(cg$membershi==tg_clusters_all[i])))
	}
}
tg_genes_all<-list()
for(i in 1:length(Cluster_list))
{
	tg_genes<-c()
	for(j in 1:length(Cluster_list[[i]]))
	{
		tg_genes1<-names(which(NMF_indi_all[,Cluster_list[[i]][j]]==1))
		eee<-cor(t(data_CORS_cancer[tg_genes1,]))
		diag(eee)<-0
		tg_genes1<-names(which(apply(eee,1,mean)>0.2))
		tg_genes<-c(tg_genes,tg_genes1)
	}
	tg_genes_all[[i]]<-unique(tg_genes)
}
All_genes_all<-c()
for(i in 1:length(tg_genes_all))
{
	All_genes_all<-c(All_genes_all,tg_genes_all[[i]])
}
All_genes_all<-unique(All_genes_all)
All_genes_all_label<-rep(0,length(All_genes_all))
names(All_genes_all_label)<-All_genes_all
Indi_all<-c()
for(i in 1:length(tg_genes_all))
{
	ccc<-All_genes_all_label
	ccc[tg_genes_all[[i]]]<-1
	Indi_all<-cbind(Indi_all,ccc)
}
colnames(Indi_all)<-nn
return(Indi_all)
}

BCV_ttest2<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
{
          x<-data0
        fff_cc<-c()
        for(kk in 1:rounds)
        {
                cv_result <- cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
                fff_cc<-rbind(fff_cc,cv_result$msep)
        }
        pp<-c()
        ddd<-apply(fff_cc,2,mean)
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
##################
merge_R1_2<-function(RBases_R1_2,R1_list,cor_cut=0.9)
{
	cor_c<-cor(t(RBases_R1_2))
	h<-hclust(dist(cor(t(RBases_R1_2))))
	N<-1
	tc<-cutree(h,N)
	st<-1
	while(st>0)
	{
		N<-N+1
		tc1<-cutree(h,N)
		tg_ids<-unique(tc1)
		st<-0
		for(j in 1:length(tg_ids))
		{
			st<-st+sum(cor_c[which(tc1==tg_ids[j]),which(tc1==tg_ids[j])]<cor_cut)
		}
	}
	cc<-tc1
	tg_ids<-unique(tc1)
	st<-0
	merged_R1_2_list<-list()
	for(j in 1:length(tg_ids))
	{
		tg_id_c<-which(tc1==tg_ids[j])
		#print(j)
		#print(tg_id_c)
		#print(names(R1_list)[tg_id_c])
		cc<-c()
		for(k in 1:length(tg_id_c))
		{
			cc<-c(cc,R1_list[[tg_id_c[k]]])
		}
		merged_R1_2_list[[j]]<-sort(unique(cc))
	}
	return(merged_R1_2_list)
}
