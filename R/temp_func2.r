Step2plus_Celltype_marker_inference_new_LowReso<-function(tg_R1_lists,data_CORS_cancer,tg_R1_cut=8,tg_R1_list_stat,cell_type_enrich_cut=0.5,resolution_level=resolution_level0,hcutn0=40)
{
        #have two sets of parameters for BCV t tests
        #1. For testing total rank and determining adding the base of a new cell type (strong condition)
        #2. For the linear model based cell type specific marker selected (weak condition)
        #3. Parameter set includes: number of BCV bins, number of rounds for t test and msep cut
	  #4. Add a new step, select top K R1 bases by using hclust, and add R1 bases by a mixed score of (1) three combination linear dependency and (2) cell type uniqueness
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
	print( dim_tt)
        #remove the leave one out step and add a hclust pre-filtering step
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

        print("Linking graph based cell type selection")
        ddd<-combn(length(tgs),3)
        comp_ids<-c()
        stat_all<-rep(0,length(tgs))
        stat_2total<-rep(0,length(tgs))
        names(stat_all)<-1:length(tgs)
        names(stat_2total)<-1:length(tgs)

        data_c<-data_CORS_cancer
	
	stat_all2<-stat_all	
	stat_2total2<-stat_2total

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
	rownames(Base_all)<-names(tgs)

	
for(i in 1:ncol(ddd))
{
	  st<-0
	  indi_c<--1
	  if(length(table(IM_low_reso_group_indi[rownames(Base_all)[ddd[,i]]]))==1)
	  {
		st<-1
		indi_c<-IM_low_reso_group_indi[rownames(Base_all)[ddd[,i]]][1]
	  }
        tg_genes<-unique(c(tgs[[ddd[1,i]]],tgs[[ddd[2,i]]],tgs[[ddd[3,i]]]))
     	  pp<-sum(BCV_ttest2(data_c[tg_genes,],rounds=20,maxrank0=5,msep_cut=0.01)<0.001)
        if(pp==2)
        {
		comp_ids<-c(comp_ids,i)
            #print(i)
           lm0<-lm(Base_all[ddd[1,i],]~Base_all[ddd[2,i],]+Base_all[ddd[3,i],]+0)
		if(st==0)
		{
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
		if(st==1)
		{
                  	if(sum(is.na(coefficients(lm0)))==0)
                        {
                              stat_2total2[ddd[,i]]<-stat_2total2[ddd[,i]]+1
                        	if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==1))
                        	{
                                stat_all2[ddd[1,i]]<-stat_all2[ddd[1,i]]+1
                        	}
                        	if((sign(coefficients(lm0))[1]==-1)&(sign(coefficients(lm0))[2]==1))
                        	{
                                stat_all2[ddd[3,i]]<-stat_all2[ddd[3,i]]+1
                        	}
                        	if((sign(coefficients(lm0))[1]==1)&(sign(coefficients(lm0))[2]==-1))
                        	{
                                stat_all2[ddd[2,i]]<-stat_all2[ddd[2,i]]+1
                        	}
                        }
		}
       }
         if(i%%100==1)
        {
                print(i)
        }
}


ooo<-stat_all2[which(stat_2total2!=0)]/stat_2total2[which(stat_2total2!=0)]
sort(ooo,decreasing=T)[24:29]
names(tg_R1_lists_selected)[as.numeric(names(sort(ooo,decreasing=T)))]


ppp<-as.numeric(names(sort(stat_all2[which(stat_2total2!=0)]/stat_2total2[which(stat_2total2!=0)],decreasing=T)))[1:20]
qqq<-setdiff(1:40,ppp)

stat_all0<-stat_all/stat_2total
scores10<-stat_all0
s1<-scores10[order(-scores10)]
ss<-rownames(Base_all)[order(-scores10)]
score_plus<-rep(0,length(ss))
names(score_plus)<-names(s1)
ss0<-unique(ss)
ss0_counts<-rep(1,length(ss0))
names(ss0_counts)<-ss0
for(i in 1:length(ss))
{
	score_plus[i]<-ss0_counts[ss[i]]
	ss0_counts[ss[i]]<-ss0_counts[ss[i]]+1
}
score_plus0<-(1/score_plus)^2
s2<-score_plus0+s1
s2<-s2[order(-s2)]


total_dim<-dim_tt
N<-0
tt_selected_genes<-c()
selected_p<-c()
selected_c<-c()
tg_d_stat<-c()
selected_id<-c()
for(i in 1:length(s2))
{
        if(N<=total_dim)
        {
                tg_id_c<-as.numeric(names(s2)[i])
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
#add_R1_base_full<-list()
if(length(tg_2badded_cells)>0)
{
        ccc<-rep(0,length(tg_2badded_cells))
        names(ccc)<-tg_2badded_cells
        ttt<-s2
        for(i in 1:length(ttt))
        {
                if(sum(names(tgs)[as.numeric(names(ttt[i]))]==tg_2badded_cells))
                {
                        tg<-names(tgs)[as.numeric(names(ttt[i]))]
                        if(ccc[tg]==0)
                        {
                                add_R1_base_cut[[which(tg_2badded_cells==tg)]]<-tgs[[as.numeric(names(ttt[i]))]]
                                #add_R1_base_full[[which(tg_2badded_cells==tg)]]<-tg_R1_lists[[as.numeric(names(ttt[i]))]]
                                ccc[tg]<-1
                        }
                }
        }
}
names(add_R1_base_cut)<-tg_2badded_cells

#############################
nn<-rep("",length(tg_R1_list_stat_selected))
for(i in 1:length(tg_R1_list_stat_selected))
{
        if(max(tg_R1_list_stat_selected[[i]][min(tg_R1_cut,nrow(tg_R1_list_stat_selected[[i]])),resolution_level])>cell_type_enrich_cut)
        {
                ccc<-tg_R1_list_stat_selected[[i]][min(tg_R1_cut,nrow(tg_R1_list_stat_selected[[i]])),resolution_level]
                nn[i]<-names(which(ccc==max(ccc))[1])
        }
}
selected_cell_types<-intersect(resolution_level,unique(nn))

cellreso_R1_base_cut<-list()
#cellreso_R1_base_full<-list()

for(i in 1:length(selected_cell_types))
{
        cc1<-c()
        #cc2<-c()
        for(j in 1:length(nn))
        {
                if(nn[j]==selected_cell_types[i])
                {
                        cc1<-c(cc1,tgs[[j]])
                        #cc2<-c(cc2,tg_R1_lists[[j]])
                }
        }
        cellreso_R1_base_cut[[i]]<-unique(cc1)
        #cellreso_R1_base_full[[i]]<-unique(cc2)
}
names(cellreso_R1_base_cut)<-selected_cell_types
#names(cellreso_R1_base_full)<-selected_cell_types

#############################
supp_results<-list(tg_R1_lists_selected,tg_R1_list_stat_selected,selected_id,tg_d_stat,dim_tt,stat_all,stat_2total,hclust_R1_screen,Base_hclust_screen_selected,R_bases_hclust_top_correlations)
names(supp_results)<-c("tg_R1_lists_selected","tg_R1_list_stat_selected","selected_id","tg_d_stat","total_cell_type_number","stat_all","stat_2total","hclust_R1_screen","Base_hclust_screen_selected","R_bases_hclust_top_correlations")

selected_R1_base_cut<-list()
#selected_R1_base_full<-list()

for(i in 1:length(selected_id))
{
        selected_R1_base_cut[[i]]<-tgs[[selected_id[i]]]
        #selected_R1_base_full[[i]]<-tg_R1_lists[[selected_id[i]]]
}
names(selected_R1_base_cut)<-names(tgs)[selected_id]
#names(selected_R1_base_full)<-names(tgs)[selected_id]

R2_selected_cell_type_markers<-list(selected_R1_base_cut,add_R1_base_cut,cellreso_R1_base_cut,supp_results)
names(R2_selected_cell_type_markers)<-c("selected_R1_base_cut","add_R1_base_cut","cellreso_R1_base_cut","supp_results")
return(R2_selected_cell_type_markers)
}
