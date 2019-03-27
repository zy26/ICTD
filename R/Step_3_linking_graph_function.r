
BCV_ttest3<-function(data0,rounds=20,slice0=2,maxrank0=30,msep_cut=0.001)
{
          x<-data0
        fff_cc<-c()
        for(kk in 1:rounds)
        {
                cv_result <- bcv::cv.svd.gabriel(x, slice0, slice0, maxrank = maxrank0)
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
        return(list(pp,fff_cc))       
}

resolution_level0<-c("fibroblast","adipocytes","neuron","endothelial","Bcell","NKcell","CD4Tcell","CD8Tcell" ,"dendritic","monocytes","macrophage","neutrophils")

Step2plus_Celltype_marker_inference<-function(tg_R1_lists,data_CORS_cancer,tg_R1_cut=6,tg_R1_list_stat,cell_type_enrich_cut=0.5,resolution_level=resolution_level0)
{
	#have two sets of parameters for BCV t tests
	#1. For testing total rank and determining adding the base of a new cell type (strong condition)
	#2. For the linear model based cell type specific marker selected (weak condition)
	#3. Parameter set includes: number of BCV bins, number of rounds for t test and msep cut
	tg_R1_lists_st<-tg_R1_lists
	for(i in 1:length(tg_R1_lists_st))
	{
		tg_R1_lists_st[[i]]<-tg_R1_lists_st[[i]][1:min(length(tg_R1_lists_st[[i]]),tg_R1_cut)]
	}
	
	print("Compute_total_rank!")
	tgs<-tg_R1_lists_st
	tg_all_genes<-c()
	for(i in 1:length(tgs))
	{
		tg_all_genes<-c(tg_all_genes,tgs[[i]])
	}
	tg_all_genes<-unique(tg_all_genes)
	tg_data_ccc<-data_CORS_cancer[tg_all_genes,]
	BCV_stat_c<-BCV_ttest3(tg_data_ccc,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
	dim_tt<-sum(BCV_stat_c[[1]]<0.01)

	print("Compute leave one base out residual rank!")	#time consuming
	BCV_stat_leave1out<-list()
	BCV_p_leave1out<-c()
	for(i in 1:length(tg_R1_lists))
	{
		print(i)
		tg_all_genes0<-c()
		for(j in 1:length(tgs))
		{
			if(j!=i)
			{
				tg_all_genes0<-c(tg_all_genes0,tgs[[j]])
			}
		}
		tg_all_genes0<-unique(tg_all_genes0)
		tg_data_ccc0<-data_CORS_cancer[tg_all_genes0,]		
		BCV_stat_c0<-BCV_ttest3(tg_data_ccc0,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
		dim_tt0<-sum(BCV_stat_c0[[1]]<0.01)
		BCV_stat_leave1out[[i]]<-BCV_stat_c0
		BCV_p_leave1out<-c(BCV_p_leave1out,dim_tt0)
	}

	print("Linking graph based cell type selection")	#time consuming
	ddd<-combn(length(tgs),3)
	comp_ids<-c()
	stat_all<-rep(0,length(tgs))
	stat_2total<-rep(0,length(tgs))
	names(stat_all)<-1:length(tgs)
	data_c<-data_CORS_cancer

	Base_all<-c()
	for(i in 1:length(tgs))
	{
		tg_data_c<-data_c[tgs[[i]],]
		cc<-svd(tg_data_c)$v[,1]
		ccc<-cor(cc,t(tg_data_c))
		if(mean(ccc)<0)
  	    {
			cc<--cc
		}
		Base_all<-rbind(Base_all,cc)
	}

	for(i in 1:ncol(ddd))
	{
     		tg_genes<-unique(c(tgs[[ddd[1,i]]],tgs[[ddd[2,i]]],tgs[[ddd[3,i]]]))
     		pp<-sum(BCV_ttest2(data_c[tg_genes,],rounds=20,maxrank0=5)<0.001)
     		if(pp==2)
    		{
				comp_ids<-c(comp_ids,i)
                        #print(i)
                lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
                if(sum(is.na(coefficients(lm0)))==0)
				{
					stat_2total[ddd[,i]]<-stat_2total[ddd[,i]]+1
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

			if(i%%100==1) print(i)
	}


	stat_all0<-stat_all/stat_2total

	total_dim<-dim_tt
	N<-0
	tt_selected_genes<-c()
	selected_p<-c()
	selected_c<-c()
	tg_d_stat<-c()
	selected_id<-c()
	for(i in 1:length(stat_all))
	{
		if(N<=total_dim)
		{
			tg_id_c<-order(-stat_all0)[i]
			print(c(tg_id_c,tg_id_c))
			if(length(selected_p)==0)
			{
				selected_p<-tgs[[tg_id_c]]
				tg_d_stat<-1
				selected_id<-c(tg_id_c)
				N<-N+1
				print(N)
				print(selected_id)
			}
			else
			{
				selected_c<-c(selected_p,tgs[[tg_id_c]])
				#selected_c<-c(selected_p,tgs[[6]])
				selected_c<-unique(selected_c)
				tg_data_ccc0<-data_CORS_cancer[selected_c,]		
				BCV_stat_c0<-BCV_ttest3(tg_data_ccc0,rounds=100,slice0=2,maxrank0=30,msep_cut=0.0001)
				tg_d_c<-sum(BCV_stat_c0[[1]]<0.001)#which(BCV_stat_c0[[1]]>0.001)[1]-1
				print(tg_d_c)
				if(tg_d_c>tg_d_stat[length(tg_d_stat)])
				{
					tg_d_stat<-c(tg_d_stat,tg_d_c)
					selected_id<-c(selected_id,tg_id_c)
					selected_p<-selected_c
					N<-N+1
					print(N)
					print(selected_id)
				}
				print("")
			} 
		}
	}


	tg_2badded_cells<-setdiff(unique(names(tgs)),unique(names(tgs)[selected_id]))
	add_R1_base_cut<-list()
	add_R1_base_full<-list()
	if(length(tg_2badded_cells)>0)
	{
		ccc<-rep(0,length(tg_2badded_cells))
		names(ccc)<-tg_2badded_cells
		ttt<-order(-stat_all0)
		for(i in 1:length(ttt))
		{
			if(sum(names(tgs)[ttt[i]]==tg_2badded_cells))
			{
				tg<-names(tgs)[ttt[i]]
				if(ccc[tg]==0)
				{
					add_R1_base_cut[[which(tg_2badded_cells==tg)]]<-tgs[[ttt[i]]]
					add_R1_base_full[[which(tg_2badded_cells==tg)]]<-tg_R1_lists[[ttt[i]]]
					ccc[tg]<-1
				}
			}
		}
	}
	names(add_R1_base_cut)<-tg_2badded_cells
	names(add_R1_base_full)<-tg_2badded_cells


	#############################
	nn<-rep("",length(tg_R1_list_stat))
	for(i in 1:length(tg_R1_list_stat))
	{
		if(max(tg_R1_list_stat[[i]][min(tg_R1_cut,nrow(tg_R1_list_stat[[i]])),resolution_level])>cell_type_enrich_cut)
		{
			ccc<-tg_R1_list_stat[[i]][min(tg_R1_cut,nrow(tg_R1_list_stat[[i]])),resolution_level]
			nn[i]<-names(which(ccc==max(ccc))[1])
		}
	}
	selected_cell_types<-intersect(resolution_level,unique(nn))

	cellreso_R1_base_cut<-list()
	cellreso_R1_base_full<-list()

	for(i in 1:length(selected_cell_types))
	{
		cc1<-c()
		cc2<-c()
		for(j in 1:length(nn))
		{
			if(nn[j]==selected_cell_types[i])
			{
				cc1<-c(cc1,tgs[[j]])
				cc2<-c(cc2,tg_R1_lists[[j]])
			}
		}
		cellreso_R1_base_cut[[i]]<-unique(cc1)
		cellreso_R1_base_full[[i]]<-unique(cc2)
	}
	names(cellreso_R1_base_cut)<-selected_cell_types
	names(cellreso_R1_base_full)<-selected_cell_types

	#############################
	supp_results<-list(selected_id,tg_d_stat,BCV_p_leave1out,BCV_stat_leave1out,dim_tt,stat_all,stat_2total)
	names(supp_results)<-c("selected_id","tg_d_stat","BCV_p_leave1out","BCV_stat_leave1out","dim_tt","stat_all","stat_2total")

	selected_R1_base_cut<-list()
	selected_R1_base_full<-list()

	for(i in 1:length(selected_id))
	{
		selected_R1_base_cut[[i]]<-tgs[[selected_id[i]]]
		selected_R1_base_full[[i]]<-tg_R1_lists[[selected_id[i]]]
	}
	names(selected_R1_base_cut)<-names(tgs)[selected_id]
	names(selected_R1_base_full)<-names(tgs)[selected_id]

	R2_selected_cell_type_markers<-list(selected_R1_base_full,selected_R1_base_cut,add_R1_base_cut,add_R1_base_full,cellreso_R1_base_cut,cellreso_R1_base_full,supp_results)
	names(R2_selected_cell_type_markers)<-c("selected_R1_base_full","selected_R1_base_cut","add_R1_base_cut","add_R1_base_full",
				"cellreso_R1_base_cut","cellreso_R1_base_full","supp_results")
	return(R2_selected_cell_type_markers)
}


################################

Building_NMF_input_new<-function(tg_data,tg_list,tg_list_add=list(),module_selected)
{
S_indi<-c(rep(1,length(tg_list)),rep(2,length(tg_list_add)),rep(0,length(module_selected[[1]])))
dd<-c(names(tg_list),names(tg_list_add))
nn<-c(paste(dd,1:length(dd),sep="_"),paste("Cancer",1:length(module_selected[[1]]),sep="_"))
names(S_indi)<-nn
Cancer_Base_all<-c()
for(i in 1:length(module_selected[[1]]))
{
	tg_data_c<-tg_data[module_selected[[1]][[i]],]
	cc<-svd(tg_data_c)$v[,1]
	ccc<-cor(cc,t(tg_data_c))
	if(mean(ccc)<0)
      {
		cc<--cc
	}
	Cancer_Base_all<-rbind(Cancer_Base_all,cc)
}
pp_pre<-matrix(0,length(S_indi),ncol(tg_data))
pp_pre[which(S_indi==0),]<-Cancer_Base_all
colnames(pp_pre)<-colnames(tg_data)
rownames(pp_pre)<-nn
tg_all_genes<-c()
for(i in 1:length(tg_list))
{
	tg_all_genes<-c(tg_all_genes,tg_list[[i]])
}
if(length(tg_list_add)>0)
{
	for(i in 1:length(tg_list_add))
	{
		tg_all_genes<-c(tg_all_genes,tg_list_add[[i]])
	}
}
tg_all_genes<-unique(tg_all_genes)
NMF_indi_all<-matrix(0,length(tg_all_genes),length(S_indi))
colnames(NMF_indi_all)<-names(S_indi)
rownames(NMF_indi_all)<-tg_all_genes
N<-0
for(i in 1:length(tg_list))
{
	N<-N+1
	NMF_indi_all[tg_list[[i]],N]<-1
}
if(length(tg_list_add)>0)
{
	for(i in 1:length(tg_list_add))
	{
		N<-N+1
		NMF_indi_all[tg_list_add[[i]],N]<-1
	}
}
NMF_data<-tg_data[rownames(NMF_indi_all),]
NMF_P_pre<-pp_pre
NMF_inputs<-list(NMF_data,NMF_indi_all,NMF_P_pre,S_indi)
names(NMF_inputs)<-c("NMF_data","NMF_indi_all","NMF_P_pre","S_indi")
return(NMF_inputs)
}

###################
Building_NMF_input_no_cancer <- function(tg_data,tg_list,tg_list_add=list()){

S_indi<-c(rep(1,length(tg_list)),rep(2,length(tg_list_add)))

dd<-c(names(tg_list),names(tg_list_add))
nn<-c(paste(dd,1:length(dd),sep="_"))
names(S_indi)<-nn

pp_pre<-matrix(0,length(S_indi),ncol(tg_data))

colnames(pp_pre)<-colnames(tg_data)
rownames(pp_pre)<-nn

tg_all_genes<-c()
for(i in 1:length(tg_list))
{
	tg_all_genes<-c(tg_all_genes,tg_list[[i]])
}
if(length(tg_list_add)>0)
{
	for(i in 1:length(tg_list_add))
	{
		tg_all_genes<-c(tg_all_genes,tg_list_add[[i]])
	}
}
tg_all_genes<-unique(tg_all_genes)
NMF_indi_all<-matrix(0,length(tg_all_genes),length(S_indi))
colnames(NMF_indi_all)<-names(S_indi)
rownames(NMF_indi_all)<-tg_all_genes

N<-0
for(i in 1:length(tg_list))
{
	N<-N+1
	NMF_indi_all[tg_list[[i]],N]<-1
}
if(length(tg_list_add)>0)
{
	for(i in 1:length(tg_list_add))
	{
		N<-N+1
		NMF_indi_all[tg_list_add[[i]],N]<-1
	}
}

NMF_data<-tg_data[rownames(NMF_indi_all),]
NMF_P_pre<-pp_pre

NMF_inputs<-list(NMF_data,NMF_indi_all,NMF_P_pre,S_indi)
names(NMF_inputs)<-c("NMF_data","NMF_indi_all","NMF_P_pre","S_indi")

return(NMF_inputs)

}
