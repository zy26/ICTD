#clean_rank1_module_new
#Step2plus_Celltype_marker_inference_fixedCT
#Step2plus_Celltype_marker_inference_Hclust
#Step2plus_Celltype_marker_inference_Hclust_plus
#Step2plus_Celltype_marker_inference_Hclust_CTES
#Explanation_BASE
#Extract_top_base_cor_ID
#Compute_Rbase_SVD





clean_rank1_module_new<-function(data_c,module_info,module_rank,st0=8,RR=50,msep_cut0=0.01)
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
                                pp<-BCV_ttest2(data_c[tg_genes,],rounds=RR,maxrank0=5,msep_cut=msep_cut0)
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

R1_list_filtering_step1_new<-function(list_c2,data_CORS_cancer,max_cut=20,cutn0=20,cut10=0.8,IM_id_list,immune_cell_uni_table=immune_cell_uni_table0_GPL570)
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

pp_all<-c()
for(i in 1:length(tg_marker_lists))
{
        pp<-sum(BCV_ttest2(data_CORS_cancer[tg_marker_lists[[i]],],maxrank0=20,msep_cut=0.01)<0.001)     
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f1<-clean_rank1_module_new(data_CORS_cancer,tg_marker_lists,pp_all,st0=6)
pp_R1_marker_list_f1.5<-cut_modules(pp_R1_marker_list_f1,cutn=cutn0)

stat_p1<-list()
for(i in 1:length(pp_R1_marker_list_f1.5))
{
        tg_gene_c<-pp_R1_marker_list_f1.5[[i]][-1]
        stat_p1[[i]]<-tg_1_rank_markers0[[names(pp_R1_marker_list_f1.5)[i]]][tg_gene_c,]
}

names(tg_marker_lists)<-names(pp_R1_marker_list_f1.5)
tg_genes_all<-c()
for(i in 1:length(pp_R1_marker_list_f1.5))
{
        tg_genes_all<-c(tg_genes_all,pp_R1_marker_list_f1.5[[i]])
}
tg_genes_all<-unique(sort(tg_genes_all))

print("filter 1 and stat done!")
#ccc<-MAP_GPL570_genes2(R1_markers_f1)
#names(ccc)<-names(Top_cell_proportion)
#table(names(Top_cell_proportion))
R1_markers_f1<-pp_R1_marker_list_f1.5

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

R1_marker_list_f2.5_stat<-rank_based_module_sorting(data_CORS_cancer,R1_marker_list_f2,IM_id_list,immune_cell_uni_table=immune_cell_uni_table)
R1_marker_list_f2.5<-R1_marker_list_f2.5_stat[[1]]

pp_all<-c()
for(i in 1:length(R1_marker_list_f2.5))
{
        pp<-sum(BCV_ttest2(data_CORS_cancer[R1_marker_list_f2.5[[i]],],maxrank0=20,msep_cut=0.01)<0.001)     
        pp_all<-c(pp_all,pp)
}
pp_R1_marker_list_f3<-clean_rank1_module(data_CORS_cancer,R1_marker_list_f2.5,pp_all,st0=6)
pp_R1_marker_list_f3.5<-cut_modules(pp_R1_marker_list_f3,cutn=cutn0)
R1_marker_list_f3.5_stat<-rank_based_module_sorting(data_CORS_cancer,pp_R1_marker_list_f3.5,IM_id_list,immune_cell_uni_table=immune_cell_uni_table)

print("filter 2 done!")
ccc<-c()
nn<-c()
for(i in 1:length(pp_R1_marker_list_f3.5))
{       
        ccc0<-c()
        for(j in 1:length(IM_id_list))
        {       
      if(length(IM_id_list[[j]])>1)
      {
            cc0<-apply(immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]],IM_id_list[[j]]],1,sum)/sum((1/(1:length(IM_id_list[[j]]))))
      }
      else
      {
           cc0<-immune_cell_uni_table[pp_R1_marker_list_f3.5[[i]],IM_id_list[[j]]]
      }
                ccc0<-cbind(ccc0,cc0)
        }
        colnames(ccc0)<-names(IM_id_list)
      ddd<-apply(ccc0,2,mean)
      ccc<-rbind(ccc,ddd)
      nn<-c(nn,colnames(ccc0)[which(ddd==max(ddd))[1]])
}
rownames(ccc)<-nn
cell_enrich_stat<-ccc
names(pp_R1_marker_list_f3.5)<-nn
rrr<-rep(1,length(pp_R1_marker_list_f3.5))
Filter_1_result_list<-list(pp_R1_marker_list_f1,R1_marker_list_f2,R1_marker_list_f3.5_stat,pp_R1_marker_list_f3.5,rrr,cell_enrich_stat)
names(Filter_1_result_list)<-c("R1_marker_list_f1","R1_marker_list_f2","R1_marker_list_f3.5_stat","R1_marker_list_f3.5","R1_marker_list_rank","R1_marker_list_f3.5_cell_enrich_stat")
return(Filter_1_result_list)
}

Step2plus_Celltype_marker_inference_fixedCT<-function(tg_R1_lists,tg_R1_list_stat)
{
	stat_c<-tg_R1_list_stat
	cc<-c()
	nn<-c()
	for(i in 1:length(stat_c))
	{
	aaa<-apply(stat_c[[i]],2,mean)
	cc<-c(cc,max(aaa))
	nn<-c(nn,names(which(aaa==max(aaa))[1]))
	}
	names(cc)<-nn
	dd<-unique(sort(nn))
	tg_id_selected<-c()
	cell_enrich_s<-c()
	for(i in 1:length(dd))
	{
		tg_id_c<-which(names(cc)==dd[i])
		tg_id_selected<-c(tg_id_selected,tg_id_c[order(-cc[tg_id_c])[1]])
		cell_enrich_s<-c(cell_enrich_s,cc[tg_id_c[order(-cc[tg_id_c])[1]]])
	}
	tg_R1_selected<-list()
	for(i in 1:length(tg_id_selected))
	{
	       tg_R1_selected[[i]]<-tg_R1_lists[[tg_id_selected[i]]]
	}
	names(tg_R1_selected)<-dd
	names(cell_enrich_s)<-dd
	cc<-list(tg_R1_selected,cell_enrich_s)
	names(cc)<-c("tg_R1_selected_fixedCT","tg_R1_selected_fixedCT")
	return(cc)
}


Step2plus_Celltype_marker_inference_Hclust<-function(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat,cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40)
{
	 tg_R1_lists_st<-tg_R1_lists
        for(i in 1:length(tg_R1_lists_st))
        {
                tg_R1_lists_st[[i]]<-tg_R1_lists_st[[i]][1:min(length(tg_R1_lists_st[[i]]),tg_R1_cut)]
        }
        
        print("Compute_total_rank!")
        tg_all_genes<-c()
        for(i in 1:length(tg_R1_lists_st))
        {
                tg_all_genes<-c(tg_all_genes,tg_R1_lists_st[[i]])
        }
        tg_all_genes<-unique(tg_all_genes)
        tg_data_ccc<-data_CORS_cancer[tg_all_genes,]
        BCV_stat_c<-BCV_ttest3(tg_data_ccc,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
        dim_tt<-sum(BCV_stat_c[[1]]<0.01)
        print("Total Cell Dim")
	print(dim_tt)
	data_c<-data_CORS_cancer
	Base_all<-c()
	for(i in 1:length(tg_R1_lists_st))
	{
        tg_data_c<-data_c[tg_R1_lists_st[[i]],]
        cc<-svd(tg_data_c)$v[,1]
        ccc<-cor(cc,t(tg_data_c))
        if(mean(ccc)<0)
      {
                cc<--cc
        }
        Base_all<-rbind(Base_all,cc)
	}
	rownames(Base_all)<-names(tg_R1_lists_st)
	bbb_all0<-Base_all
	rownames(bbb_all0)<-1:nrow(bbb_all0)
	hcutn0<-min(hcutn0,nrow(bbb_all0))
	hcutn0<-min(hcutn0,dim_tt)
	data_genes_all<-data_CORS_cancer[tg_all_genes,]
	hcutn1<-min(hcutn0+3,nrow(bbb_all0))
	hclust_R1_screen<-hclust_screen_top_bases(bbb_all0,data_genes_all,hcutn=hcutn1)
	ccc<-c()
	for(i in 1:length(hclust_R1_screen[[1]]))
	{
		aaa<-cor(t(hclust_R1_screen[[1]][[i]]))
		diag(aaa)<-0
		ccc<-c(ccc,max(aaa))
	}
	R_bases_hclust_top_correlations<-ccc
	hh<-hclust_R1_screen[[3]][[hcutn0]]
	Base_hclust_screen_result<-Base_screen(bbb_all0,hh,tg_R1_list_stat)
	Base_hclust_screen_selected<-as.numeric(Base_hclust_screen_result[[1]])

	
	tg_R1_lists_selected<-list()
	tg_R1_list_stat_selected<-list()
	nn<-c()
	for(i in 1:length(Base_hclust_screen_selected))
	{
		tg_R1_lists_selected[[i]]<- tg_R1_lists_st[[Base_hclust_screen_selected[[i]]]]
		tg_R1_list_stat_selected[[i]]<-tg_R1_list_stat[[Base_hclust_screen_selected[[i]]]][tg_R1_lists_selected[[i]],]
		cc<-apply(tg_R1_list_stat_selected[[i]],2,mean)
		nn<-c(nn,names(which(cc==max(cc))[1]))
	}
	names(tg_R1_lists_selected)<-nn
	tgs<-tg_R1_lists_selected
	cc<-list(tgs,R_bases_hclust_top_correlations,dim_tt,hcutn0,Base_hclust_screen_result)
	names(cc)<-list("tg_R1_selected_Hclust","R1Base_correlation_vs_hctn","Inferred_CT_dimension","Forced_CT_number","Base_hclust_screen_result")
	return(cc)
}

Step2plus_Celltype_marker_inference_Hclust_plus<-function(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat,cell_type_enrich_cut=0.5,
resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.95,IM_reso_level=IM_reso_level)
{
        stat_c<-tg_R1_list_stat
        cc<-c()
        nn<-c()
        for(i in 1:length(stat_c))
        {
        aaa<-apply(stat_c[[i]],2,mean)
        cc<-c(cc,max(aaa))
        nn<-c(nn,names(which(aaa==max(aaa))[1]))
        }
        names(cc)<-nn
        dd<-unique(sort(nn))
        tg_id_selected<-c()
        cell_enrich_s<-c()
        for(i in 1:length(dd))
        {
                tg_id_c<-which(names(cc)==dd[i])
                tg_id_selected<-c(tg_id_selected,tg_id_c[order(-cc[tg_id_c])[1]])
                cell_enrich_s<-c(cell_enrich_s,cc[tg_id_c[order(-cc[tg_id_c])[1]]])
        }
        tg_R1_selected<-list()
        for(i in 1:length(tg_id_selected))
        {
               tg_R1_selected[[i]]<-tg_R1_lists[[tg_id_selected[i]]]
        }
        names(tg_R1_selected)<-dd
        names(cell_enrich_s)<-dd
        dd0<-dd
        names(dd0)<-tg_id_selected      

         tg_R1_lists_st<-tg_R1_lists
       for(i in 1:length(tg_R1_lists_st))
       {
            tg_R1_lists_st[[i]]<-tg_R1_lists_st[[i]][1:min(length(tg_R1_lists_st[[i]]),tg_R1_cut)]
       }

        print("Compute_total_rank!")
        tg_all_genes<-c()
        for(i in 1:length(tg_R1_lists_st))
        {
                tg_all_genes<-c(tg_all_genes,tg_R1_lists_st[[i]])
        }
        tg_all_genes<-unique(tg_all_genes)
        tg_data_ccc<-data_CORS_cancer[tg_all_genes,]
        BCV_stat_c<-BCV_ttest3(tg_data_ccc,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
        dim_tt<-sum(BCV_stat_c[[1]]<0.01)
        print("Total Cell Dim")
        print(dim_tt)
        data_c<-data_CORS_cancer
        Base_all<-c()
        for(i in 1:length(tg_R1_lists_st))
        {
        tg_data_c<-data_c[tg_R1_lists_st[[i]],]
        cc<-svd(tg_data_c)$v[,1]
        ccc<-cor(cc,t(tg_data_c))
        if(mean(ccc)<0)
        {
                cc<--cc
        }
        Base_all<-rbind(Base_all,cc)
        }
        rownames(Base_all)<-names(tg_R1_lists_st)
        bbb_all0<-Base_all
        rownames(bbb_all0)<-1:nrow(bbb_all0)
        hcutn0<-min(hcutn0,nrow(bbb_all0))
        hcutn0<-min(hcutn0,dim_tt)
        data_genes_all<-data_CORS_cancer[tg_all_genes,]
        hcutn1<-min(hcutn0+3,nrow(bbb_all0))
        hclust_R1_screen<-hclust_screen_top_bases(bbb_all0,data_genes_all,hcutn=hcutn1)
       
	 R1_hhh<-hclust_R1_screen[[3]]
        ccc<-c()
        for(i in 1:length(hclust_R1_screen[[1]]))
        {
                aaa<-cor(t(hclust_R1_screen[[1]][[i]]))
                diag(aaa)<-0
                ccc<-c(ccc,max(aaa))
        }
        R_bases_hclust_top_correlations<-ccc
        hcutn2<-max(which(R_bases_hclust_top_correlations<hclust_cor_cut))
        hh<-hclust_R1_screen[[3]][[hcutn2]]
        Base_hclust_screen_result<-Base_screen(bbb_all0,hh,tg_R1_list_stat)
        Base_hclust_screen_selected<-as.numeric(Base_hclust_screen_result[[1]])
        ttt<-unique(hh)
        Base_hclust_screen_selected_hp<-Base_hclust_screen_selected
        IM_reso_level0<-IM_reso_level

	 for(i in 1:length(ttt))
        {
                cc1<-c()
                cc1<-intersect(names(which(hh==ttt[i])),Base_hclust_screen_result[[1]][i])
                cc2<-c()
                cc2<-intersect(as.character(tg_id_selected),names(which(hh==ttt[i])))
                if(length(cc2)==1)
                {
                        Base_hclust_screen_selected_hp[i]<-as.numeric(cc2)
                }
                if(length(cc2)>1)
                {
                        ccc<-IM_reso_level[dd0[cc2]]
                        cc3<-cc2[which(ccc==max(ccc))]
                        Base_hclust_screen_selected_hp[i]<-as.numeric(cc3)[1]
                }
        }
        print("done!")

       
        
        tg_R1_lists_selected<-list()
        tg_R1_list_stat_selected<-list()
        nn<-c()
        for(i in 1:length(Base_hclust_screen_selected_hp))
        {
                tg_R1_lists_selected[[i]]<- tg_R1_lists_st[[Base_hclust_screen_selected_hp[[i]]]]
                tg_R1_list_stat_selected[[i]]<-tg_R1_list_stat[[Base_hclust_screen_selected_hp[[i]]]][tg_R1_lists_selected[[i]],]
                cc<-apply(tg_R1_list_stat_selected[[i]],2,mean)
                nn<-c(nn,names(which(cc==max(cc))[1]))
        }
        names(tg_R1_lists_selected)<-nn
        tgs<-tg_R1_lists_selected
        cc<-list(tgs,tg_R1_list_stat_selected,dim_tt,length(Base_hclust_screen_selected_hp),Base_hclust_screen_result)
        names(cc)<-list("tg_R1_selected_Hclust","tg_R1_list_stat_selected","Inferred_CT_dimension","Forced_CT_number","Base_hclust_screen_result")
        return(cc)
}	

###########################

Step2plus_Celltype_marker_inference_Hclust_CTES<-function(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat,cell_type_enrich_cut=0.5,
resolution_level=resolution_level0,hcutn0=40,hclust_cor_cut=0.95,IM_reso_level=IM_reso_level,extra_ctn=3)
{
        stat_c<-tg_R1_list_stat
        cc<-c()
        nn<-c()
        for(i in 1:length(stat_c))
        {
        aaa<-apply(stat_c[[i]],2,mean)
        cc<-c(cc,max(aaa))
        nn<-c(nn,names(which(aaa==max(aaa))[1]))
        }
        names(cc)<-nn
        dd<-unique(sort(nn))
        tg_id_selected<-c()
        cell_enrich_s<-c()
        for(i in 1:length(dd))
        {
                tg_id_c<-which(names(cc)==dd[i])
                tg_id_selected<-c(tg_id_selected,tg_id_c[order(-cc[tg_id_c])[1]])
                cell_enrich_s<-c(cell_enrich_s,cc[tg_id_c[order(-cc[tg_id_c])[1]]])
        }
        tg_R1_selected<-list()
        for(i in 1:length(tg_id_selected))
        {
               tg_R1_selected[[i]]<-tg_R1_lists[[tg_id_selected[i]]]
        }
        names(tg_R1_selected)<-dd
        names(cell_enrich_s)<-dd
        dd0<-dd
        names(dd0)<-tg_id_selected

       tg_R1_lists_st<-tg_R1_lists
       for(i in 1:length(tg_R1_lists_st))
       {
            tg_R1_lists_st[[i]]<-tg_R1_lists_st[[i]][1:min(length(tg_R1_lists_st[[i]]),tg_R1_cut)]
       }

       # print("Compute_total_rank!")
        tg_all_genes<-c()
        for(i in 1:length(tg_R1_lists_st))
        {
                tg_all_genes<-c(tg_all_genes,tg_R1_lists_st[[i]])
        }
        tg_all_genes<-unique(tg_all_genes)
        tg_data_ccc<-data_CORS_cancer[tg_all_genes,]
        BCV_stat_c<-BCV_ttest3(tg_data_ccc,rounds=100,slice0=4,maxrank0=30,msep_cut=0.0001)
        dim_tt<-sum(BCV_stat_c[[1]]<0.01)
        
       # print("Total Cell Dim")
       # print(dim_tt)

	  NN_ES_table_R4<-Explanation_BASE(data_CORS_cancer,tg_R1_lists_st)
	  tg_ES_scores<-diag(NN_ES_table_R4)
	
        data_c<-data_CORS_cancer
        Base_all<-c()
        for(i in 1:length(tg_R1_lists_st))
        {
        tg_data_c<-data_c[tg_R1_lists_st[[i]],]
        cc<-svd(tg_data_c)$v[,1]
        ccc<-cor(cc,t(tg_data_c))
        if(mean(ccc)<0)
        {
                cc<--cc
        }
        Base_all<-rbind(Base_all,cc)
        }
        rownames(Base_all)<-names(tg_R1_lists_st)
        bbb_all0<-Base_all
        rownames(bbb_all0)<-1:nrow(bbb_all0)
        hcutn0<-min(hcutn0,nrow(bbb_all0))
        hcutn0<-min(hcutn0,dim_tt)
        data_genes_all<-data_CORS_cancer[tg_all_genes,]
        hcutn1<-min(hcutn0+extra_ctn,nrow(bbb_all0))
        hclust_R1_screen<-hclust_screen_top_bases(bbb_all0,data_genes_all,hcutn=hcutn1)
       
	 R1_hhh<-hclust_R1_screen[[3]]
        ccc<-c()
        for(i in 1:length(hclust_R1_screen[[1]]))
        {
                aaa<-cor(t(hclust_R1_screen[[1]][[i]]))
                diag(aaa)<-0
                ccc<-c(ccc,max(aaa))
        }
        R_bases_hclust_top_correlations<-ccc
        hcutn2<-max(which(R_bases_hclust_top_correlations<hclust_cor_cut))
        hh<-hclust_R1_screen[[3]][[hcutn2]]

        Base_hclust_screen_result<-Base_screen(bbb_all0,hh,tg_R1_list_stat)
        tg_selected_Bases<-c()
	  for(i in 1:length(Base_hclust_screen_result[[2]]))
	  {
		tg_ids_c<-as.numeric(Base_hclust_screen_result[[2]][[i]][[1]])
		tg_selected_Bases<-c(tg_selected_Bases,tg_ids_c[which(tg_ES_scores[tg_ids_c]==min(tg_ES_scores[tg_ids_c]))[1]])
        }
	  Base_screen_selected<- tg_selected_Bases
        
	  selected_cell_types<-intersect(names(IM_reso_level),unique(nn))
	  rest_cell_types<-setdiff(selected_cell_types,unique(nn[Base_screen_selected]))
	  names(Base_screen_selected)<-nn[Base_screen_selected]
	  ccc<-c()
	 if(length(rest_cell_types)>0)
	{
	 for(i in 1:length(rest_cell_types))
	 {
        	tg_ids_c<-which(nn==rest_cell_types[i])
       	ccc<-c(ccc,tg_ids_c[which(tg_ES_scores[tg_ids_c]==min(tg_ES_scores[tg_ids_c]))[1]])
	 }
	 names(ccc)<-rest_cell_types
	}
       
	R1_list_c1<-list()
	for(i in 1:length(Base_screen_selected))
	{
		R1_list_c1[[i]]<-tg_R1_lists_st[[Base_screen_selected[i]]]
	}
	names(R1_list_c1)<-names(Base_screen_selected)

	R1_list_c2<-list()
	if(length(ccc)>0)
	{
	for(i in 1:length(ccc))
	{
		R1_list_c2[[i]]<-tg_R1_lists_st[[ccc[i]]]
	}
	names(R1_list_c2)<-names(ccc)
	}
        cc<-list(R1_list_c1,R1_list_c2,Base_screen_selected,ccc,dim_tt,length(Base_screen_selected),Base_hclust_screen_result)
        names(cc)<-list("tg_R1_selected_ES","tg_R1_selected_extra_CT","tg_R1_selected_ES_ids","tg_R1_selected_extra_CT_ids","Inferred_CT_dimension","Forced_CT_number","Base_hclust_screen_result")
        return(cc)
}	


Explanation_BASE<-function(data_c0,tg_R1_list_c)
{
tg_R1_lists_selected<-tg_R1_list_c
tg_R1_lists_st_ccc<-tg_R1_lists_selected
data_c<-data_c0
Base_all<-c()
for(i in 1:length(tg_R1_lists_st_ccc))
{
       tg_data_c<-data_c[tg_R1_lists_st_ccc[[i]],]
       cc<-svd(tg_data_c)$v[,1]
       ccc<-cor(cc,t(tg_data_c))
       if(mean(ccc)<0)
       {
                cc<--cc
       }
       Base_all<-rbind(Base_all,cc)
}
rownames(Base_all)<-1:nrow(Base_all)
if(length(names(tg_R1_lists_selected))>1)
{
	rownames(Base_all)<-names(tg_R1_lists_selected)
}
	NN_ES_table<-matrix(0,length(tg_R1_list_c),length(tg_R1_list_c))
	for(i in 1:length(tg_R1_list_c))
	{
		tg_data_c<-data_c[tg_R1_list_c[[i]],]
		for(j in 1:length(tg_R1_list_c))
		{
			if(cor(Base_all[i,],Base_all[j,])>0)
			{
				 ttt_ccc<-Base_all[j,]%*%t(Base_all[j,])/sum((Base_all[j,])^2)
                         sss<-tg_data_c-tg_data_c%*%ttt_ccc
                         NN_ES_table[i,j]<-mean(apply(sss^2,1,mean))
			}
			else
			{
				NN_ES_table[i,j]<-1
			}
		}
	}
	rownames(NN_ES_table)<-names(tg_R1_list_c)
	colnames(NN_ES_table)<-names(tg_R1_list_c)
	return(NN_ES_table)
}

Extract_top_base_cor_ID<-function(ddd,K=5)
{
	ccc1<-c()
	ccc2<-c()
	for(i in 1:ncol(ddd))
	{
		ccc1<-cbind(ccc1,order(-ddd[,i])[1:K])
		ccc2<-cbind(ccc2,ddd[order(-ddd[,i])[1:K],i])
	}
	colnames(ccc1)<-colnames(ddd)
	colnames(ccc2)<-colnames(ddd)
	return(list(ccc1,ccc2))
}

Compute_Rbase_SVD<-function(bulk_data,tg_R1_lists_selected)
{
tg_R1_lists_st_ccc<-tg_R1_lists_selected
data_c<-bulk_data
Base_all<-c()
for(i in 1:length(tg_R1_lists_st_ccc))
{
       tg_data_c<-data_c[tg_R1_lists_st_ccc[[i]],]
       cc<-svd(tg_data_c)$v[,1]
       ccc<-cor(cc,t(tg_data_c))
       if(mean(ccc)<0)
       {
                cc<--cc
       }
       Base_all<-rbind(Base_all,cc)
}
rownames(Base_all)<-1:nrow(Base_all)
if(length(names(tg_R1_lists_selected))>1)
{
	rownames(Base_all)<-names(tg_R1_lists_selected)
}
return(Base_all)
}



