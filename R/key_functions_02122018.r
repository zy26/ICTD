#Functions
###########
#RMSE_row<-function(data0)
#log_new<-function(data_c)
#consist_NMF_markers<-function(ccc,tg_markers_c)
#consist_NMF_markers2<-function(ccc,ccc2,tg_markers_c)
#consist_NMF_markers3<-function(ccc,ccc2,tg_markers_c)
#normalize_data<-function(data0)
#normalize_data2<-function(data0)
#normalize_data3<-function(data0)
#compute_O_Rowspace<-function(data0,k)
#compute_CO_Rowspace<-function(data0,k)
#coexpression_marker_selection<-function(tg_markers,cor_Cut=0.8,cor_C0)
#plot_BCs<-function(tg_selected_BCs,qubic_input_all,qubic_results_all,tg_selected_data_list2)
#compute_plot_BCs<-function(tg_selected_data_list2,c0=1,o0=3000,f0=0.1,BC_nrow_cut0=10)
#Compute_Jaccard_mm<-function(M0,M1)
#coexpression_marker_selection2<-function(cor_C0,marker_stats1,cor_Cut_list=seq(50,90,by=2)/100)
#plot_cor_gene_lists<-function(data0,tg_list)
#analyze_QUBIC_genes_stats<-function(BC_im_depletion_gene_stats)
#compute_Jaccard<-function(data_list)
#compute_Intersect<-function(data_list)
#norm_transform<-function(data0)
#test_cell_list<-function(BC_gene_list,tg_cell_list,all_genes)
#coexpression_marker_selection<-function(tg_markers,cor_Cut=0.8,cor_C0)

##################
#load("IM_markers_02082018.RData")

####################################
RMSE_row<-function(data0)
{
	return(apply(data0^2,1,mean))
}

log_new<-function(data_c)
{
	data_c0<-data_c
	for(i in 1:nrow(data_c0))
	{
		aaa<-log(data_c[i,])
		m<-min(data_c[i,which(data_c[i,]!=0)])
		aaa[which(data_c[i,]==0)]<-log(m*0.5)
		data_c0[i,]<-aaa
	}
	return(data_c0)
}

consist_NMF_markers<-function(ccc,tg_markers_c)
{
ddd<-c()
for(i in 1:nrow(ccc))
{
ddd<-c(ddd,which(ccc[i,]==max(ccc[i,])))
}
names(ddd)<-rownames(ccc)
aaa<-matrix(0,ncol(ccc),length(tg_markers_c))
for(i in 1:ncol(ccc))
{
for(j in 1:length(tg_markers_c))
{
aaa[i,j]<-length(intersect(names(ddd)[which(ddd==i)],tg_markers_c[[j]]))
}
}
rownames(aaa)<-1:ncol(ccc)
colnames(aaa)<-names(tg_markers_c)
fff<-list()
for(i in 1:ncol(ccc))
{
fff[[i]]<-names(ddd)[which(ddd==i)]
}
names(fff)<-1:ncol(ccc)
return(list(aaa,fff))
}


consist_NMF_markers2<-function(ccc,ccc2,tg_markers_c)
{
ccc3<-apply(ccc2,1,mean)
for(i in 1:ncol(ccc))
{
	ccc[,i]<-ccc[,i]*ccc3[i]
}
ddd<-c()
for(i in 1:nrow(ccc))
{
ddd<-c(ddd,which(ccc[i,]==max(ccc[i,])))
}
names(ddd)<-rownames(ccc)
aaa<-matrix(0,ncol(ccc),length(tg_markers_c))
for(i in 1:ncol(ccc))
{
for(j in 1:length(tg_markers_c))
{
aaa[i,j]<-length(intersect(names(ddd)[which(ddd==i)],tg_markers_c[[j]]))
}
}
rownames(aaa)<-1:ncol(ccc)
colnames(aaa)<-names(tg_markers_c)
fff<-list()
for(i in 1:ncol(ccc))
{
fff[[i]]<-names(ddd)[which(ddd==i)]
}
names(fff)<-1:ncol(ccc)
return(list(aaa,fff))
}


consist_NMF_markers3<-function(ccc,ccc2,tg_markers_c)
{
ccc3<-apply(ccc2,1,median)
for(i in 1:ncol(ccc))
{
	ccc[,i]<-ccc[,i]*ccc3[i]
}
ddd<-c()
for(i in 1:nrow(ccc))
{
ddd<-c(ddd,which(ccc[i,]==max(ccc[i,])))
}
names(ddd)<-rownames(ccc)
aaa<-matrix(0,ncol(ccc),length(tg_markers_c))
for(i in 1:ncol(ccc))
{
for(j in 1:length(tg_markers_c))
{
aaa[i,j]<-length(intersect(names(ddd)[which(ddd==i)],tg_markers_c[[j]]))
}
}
rownames(aaa)<-1:ncol(ccc)
colnames(aaa)<-names(tg_markers_c)
fff<-list()
for(i in 1:ncol(ccc))
{
fff[[i]]<-names(ddd)[which(ddd==i)]
}
names(fff)<-1:ncol(ccc)
return(list(aaa,fff))
}


normalize_data<-function(data0)
{
	data1<-data0
	for(i in 1:nrow(data1))
	{
		data1[i,]<-data1[i,]/max(data1[i,])
	}
	return(data1)
}

normalize_data2<-function(data0)
{
	data1<-data0
	for(i in 1:nrow(data1))
	{
		data1[i,]<-(data1[i,]-mean(data1[i,]))/sd(data1[i,])
	}
	return(data1)
}


normalize_data3<-function(data0)
{
	data1<-data0
	for(i in 1:nrow(data1))
	{
		data1[i,]<-(data1[i,]-min(data1[i,]))
		data1[i,]<-data1[i,]/max(data1[i,])
	}
	return(data1)
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

compute_CO_Rowspace<-function(data0,k)
{
	SVD_ccc<-svd(normalize_data2(data0))
	ttt_ccc<-t(SVD_ccc$v)[1,]%*%t(t(SVD_ccc$v)[1,])/sum((t(SVD_ccc$v)[1,])^2)
	t_ccc0<-ttt_ccc*0
	for(i in (k+1):nrow(SVD_ccc$v))
	{
		t_ccc0<-t_ccc0+t(SVD_ccc$v)[i,]%*%t(t(SVD_ccc$v)[i,])/sum((t(SVD_ccc$v)[i,])^2)
	}
	return(t_ccc0)
}


coexpression_marker_selection<-function(tg_markers,cor_Cut=0.8,cor_C0)
{
	tg_B_markers<-tg_markers
	cn<-c()
N<-0
for(i in 1:length(tg_B_markers))
{
	ccc<-names(which(cor_C0[tg_B_markers[[i]],]>0.8))
	if(length(ccc)>4)
	{
		N<-N+1
		cn<-c(cn,tg_B_markers[[i]])
		tg_B_markers_cor_lists[[N]]<-ccc
	}
}
names(tg_B_markers_cor_lists)<-cn

ff_all<-c()
tg_lists<-tg_B_markers_cor_lists
for(i in 1:length(tg_lists))
{
	ccc<-extract_data_symbol(marker_stats1,tg_lists[[i]])
	aaa<-unique(rownames(ccc))
	ddd<-c()
	for(j in 1:length(aaa))
	{
		if(sum(rownames(ccc)==aaa[j])==1)
		{
			bbb<-ccc[aaa[j],]
			bbb[which(bbb==0)]<-100
			bbb<-1/bbb*(bbb<=3)
			ddd<-rbind(ddd,bbb)
		}
		else
		{
			bbb<-ccc[which(rownames(ccc)==aaa[j]),]
			bbb[which(bbb==0)]<-100
			bbb<-1/bbb*(bbb<=3)
			ddd<-rbind(ddd,apply(bbb,2,mean))
		}
	}
	rownames(ddd)<-aaa
	fff<-apply(ddd,2,mean)
	ff_all<-rbind(ff_all,fff)
}
rownames(ff_all)<-names(tg_lists)

#library(gplots)
colors = c(0:100)/100
my_palette <- grDevices::colorRampPalette(c("white", "midnightblue"))(n =100)

heatmap.2(ff_all,Rowv=F,Colv =F,scale="none",main="",
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)
	return(list(ff_all,tg_B_markers_cor_lists))
}

##################

plot_BCs<-function(tg_selected_BCs,qubic_input_all,qubic_results_all,tg_selected_data_list2){
for(ii in 1:length(tg_selected_BCs))
{
	print(tg_cancer_keys[ii])
	tg_all_genes<-c()
	for(i in 1:length(tg_selected_BCs[[ii]][[1]]))
	{
		tg_all_genes<-c(tg_all_genes,tg_selected_BCs[[ii]][[1]][[i]])
	}
	tg_all_genes<-unique(tg_all_genes)
	aaa000<-matrix(0,length(tg_selected_BCs[[ii]][[1]]),length(tg_selected_BCs[[ii]][[1]]))
	for(i in 1:length(tg_selected_BCs[[ii]][[1]]))
	{
		for(j in 1:length(tg_selected_BCs[[ii]][[1]]))
		{
			aaa000[i,j]<-1-Compute_Jaccard_mm(qubic_input_all[[ii]][tg_all_genes,tg_selected_BCs[[ii]][[2]][[i]]],qubic_input_all[[ii]][tg_all_genes,tg_selected_BCs[[ii]][[2]][[j]]])
		}
	}
	bbb000<-hclust(as.dist(aaa000))
	plot(bbb000)
	aaa1<-bbb000$order
	hmcols2 <- grDevices::colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(length(aaa1)) 
	tg_genes_all<-c()
	tg_samples_all<-c()
	tg_col_c3<-c()
	tg_row_c3<-c()
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[i]]
		tg_samples<-tg_selected_BCs[[ii]][[2]][[i]]
		tg_samples_c<-setdiff(tg_samples,tg_samples_all)
		tg_genes_c<-setdiff(tg_genes,tg_genes_all)
		tg_genes_all<-c(tg_genes_all,tg_genes_c)
		tg_samples_all<-c(tg_samples_all,tg_samples_c)
		tg_col_c3<-c(tg_col_c3,rep(hmcols2[i],length(tg_samples_c)))
		tg_row_c3<-c(tg_row_c3,rep(hmcols2[i],length(tg_genes_c)))
	}
	#library(gplots)
	colors = c(-100:100)/100
	my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n = 200)
	colors2 = c(-100:100)/100
	my_palette2 <- grDevices::colorRampPalette(c("white","white", "midnightblue"))(n = 200)
	tg_bbb0<-qubic_input_all[[ii]][tg_genes_all,tg_samples_all]
	h<-heatmap.2(tg_bbb0,Rowv=F,Colv =F,scale="none",main=tg_cancer_keys[ii],
deColors=as.character(tg_col_c3),RowSideColors=as.character(tg_row_c3),
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_bbb1<-tg_selected_data_list2[[ii]][tg_genes_all,tg_samples_all]
	h<-heatmap.2(tg_bbb1,Rowv=F,Colv =F,scale="row",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c3),RowSideColors=as.character(tg_row_c3),
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_r_genes<-setdiff(rownames(qubic_input_all[[ii]]),tg_genes_all)
	tg_r_samples<-setdiff(colnames(qubic_input_all[[ii]]),tg_samples_all)
	cor_r_genes<-c()
	cor_r_samples<-c()
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[i]]
		tg_samples<-tg_selected_BCs[[ii]][[2]][[i]]
		cor_r_genes<-cbind(cor_r_genes,apply(cor(t(tg_selected_data_list2[[ii]][tg_r_genes,tg_samples_all]),t(tg_selected_data_list2[[ii]][tg_genes,tg_samples_all])),1,mean))
		cor_r_samples<-cbind(cor_r_samples,apply(cor(tg_selected_data_list2[[ii]][tg_genes_all,tg_r_samples],tg_selected_data_list2[[ii]][tg_genes_all,tg_samples]),1,mean))
	}
	rownames(cor_r_samples)<-tg_r_samples
	rownames(cor_r_genes)<-tg_r_genes
	cor_r_samples_indi<-c()
	tg_col_c4<-c()
	tg_row_c4<-c()
	for(i in 1:nrow(cor_r_samples))
	{
		cor_r_samples_indi<-c(cor_r_samples_indi,which(cor_r_samples[i,]==max(cor_r_samples[i,]))[1])
	}
	tg_samples_all_total<-c()
	for(i in 1:length(aaa1))
	{
		tg_samples<-tg_selected_BCs[[ii]][[2]][[i]]
		tg_samples_c<-setdiff(tg_samples,tg_samples_all_total)
		tg_samples_all_total<-c(tg_samples_all_total,tg_samples_c)
		tg_col_c4<-c(tg_col_c4,rep(hmcols2[i],length(tg_samples_c)))
		tg_samples_all_total<-c(tg_samples_all_total,names(cor_r_samples[which(cor_r_samples_indi==i),i])[order(-cor_r_samples[which(cor_r_samples_indi==i),i])])
		tg_col_c4<-c(tg_col_c4,rep("white",length(names(cor_r_samples[which(cor_r_samples_indi==i),i])[order(-cor_r_samples[which(cor_r_samples_indi==i),i])])))
	}
	cor_r_genes_indi<-c()
	for(i in 1:nrow(cor_r_genes))
	{
		cor_r_genes_indi<-c(cor_r_genes_indi,which(cor_r_genes[i,]==max(cor_r_genes[i,]))[1])
	}
	tg_genes_all_total<-c()
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[i]]
		tg_genes_c<-setdiff(tg_genes,tg_genes_all_total)
		tg_genes_all_total<-c(tg_genes_all_total,tg_genes_c)
		tg_row_c4<-c(tg_row_c4,rep(hmcols2[i],length(tg_genes_c)))
		tg_genes_all_total<-c(tg_genes_all_total,names(cor_r_genes[which(cor_r_genes_indi==i),i])[order(-cor_r_genes[which(cor_r_genes_indi==i),i])])
		tg_row_c4<-c(tg_row_c4,rep("white",length(names(cor_r_genes[which(cor_r_genes_indi==i),i])[order(-cor_r_genes[which(cor_r_genes_indi==i),i])])))
	}
	tg_bbb0<-qubic_input_all[[ii]][tg_genes_all_total,tg_samples_all_total]
	h<-heatmap.2(tg_bbb0,Rowv=F,Colv =F,scale="none",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c4),RowSideColors=as.character(tg_row_c4),
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_bbb1<-tg_selected_data_list2[[ii]][tg_genes_all_total,tg_samples_all_total]
	h<-heatmap.2(tg_bbb1,Rowv=F,Colv =F,scale="row",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c4),RowSideColors=as.character(tg_row_c4),
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
}
}



compute_plot_BCs<-function(tg_selected_data_list2,c0=1,o0=3000,f0=0.1,BC_nrow_cut0=10){
qubic_results_all<-list()
qubic_input_all<-list()
for(i in 1:length(tg_selected_data_list2))
{
	print(i)
	tg_data_cc00<-tg_selected_data_list2[[i]]
	tg_data_cc00<-tg_data_cc00[which(apply(tg_data_cc00==0,1,sum)/ncol(tg_data_cc00)<0.5),]
	tg_data_cc01<-tg_data_cc00*0
	for(j in 1:nrow(tg_data_cc00))
	{
		tg_data_cc01[j,]<-(tg_data_cc00[j,]<quantile(tg_data_cc00[j,],0.3))
	}
	qubic_input_all[[i]]<-tg_data_cc01
	aaa<-biclust(tg_data_cc01, method=BCQUD(),c = c0, o = o0, f = f0)
	qubic_results_all[[i]]<-aaa
}
names(qubic_results_all)<-tg_cancer_keys
names(qubic_input_all)<-tg_cancer_keys
print("BC_done!")
tg_selected_BCs<-list()
for(ii in 1:length(tg_selected_data_list2))
{
	print(names(qubic_results_all)[ii])
	aaa<-qubic_results_all[[ii]]
	tg_data_c<-qubic_input_all[[ii]]
	tg_selected_BCs_r<-list()
	tg_selected_BCs_c<-list()
	N<-0
	for(i in 1:ncol(attributes(aaa)$RowxNumber))
	{
		tg_r_ids<-which(attributes(aaa)$RowxNumber[,i]==1)
		tg_c_ids<-which(attributes(aaa)$NumberxCol[i,]==1)
		if(length(tg_r_ids)>BC_nrow_cut0)
		{
			N<-N+1
			tg_selected_BCs_r[[N]]<-rownames(tg_data_c)[tg_r_ids]
			tg_selected_BCs_c[[N]]<-colnames(tg_data_c)[tg_c_ids]
		}
	}
	tg_selected_BCs[[ii]]<-list(tg_selected_BCs_r,tg_selected_BCs_c)
}
names(tg_selected_BCs)<-tg_cancer_keys
tg_selected_BCs_plot<-list()
tg_selected_BCs_other_sample_genes<-list()
print("BC_Select_done!\nStart_plotting\n")
for(ii in 1:length(tg_selected_BCs))
{
	tg_selected_BCs_plot_c<-list()
	tg_selected_BCs_other_sample_genes_c<-list()
	print(tg_cancer_keys[ii])
	tg_all_genes<-c()
	for(i in 1:length(tg_selected_BCs[[ii]][[1]]))
	{
		tg_all_genes<-c(tg_all_genes,tg_selected_BCs[[ii]][[1]][[i]])
	}
	tg_all_genes<-unique(tg_all_genes)
	aaa000<-matrix(0,length(tg_selected_BCs[[ii]][[1]]),length(tg_selected_BCs[[ii]][[1]]))
	for(i in 1:length(tg_selected_BCs[[ii]][[1]]))
	{
		for(j in 1:length(tg_selected_BCs[[ii]][[1]]))
		{
			aaa000[i,j]<-1-Compute_Jaccard_mm(qubic_input_all[[ii]][tg_all_genes,tg_selected_BCs[[ii]][[2]][[i]]],qubic_input_all[[ii]][tg_all_genes,tg_selected_BCs[[ii]][[2]][[j]]])
		}
	}
	aaa1<-1:length(tg_selected_BCs[[ii]][[1]])
	if(length(tg_selected_BCs[[ii]][[1]])>2)
	{
		bbb000<-hclust(as.dist(aaa000))
		plot(bbb000)
		aaa1<-bbb000$order
	}
	hmcols2 <- grDevices::colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(length(aaa1)) 
	tg_genes_all<-c()
	tg_samples_all<-c()
	tg_col_c3<-c()
	tg_row_c3<-c()
	tg_selected_BCs_plot_gene_c<-list()
	tg_selected_BCs_plot_sample_c<-list()
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[aaa1[i]]]
		tg_samples<-tg_selected_BCs[[ii]][[2]][[aaa1[i]]]
		tg_samples_c<-setdiff(tg_samples,tg_samples_all)
		tg_genes_c<-setdiff(tg_genes,tg_genes_all)
		if(length(tg_genes_c)>0)
		{
			tg_selected_BCs_plot_gene_c[[i]]<-tg_genes_c
		}
		else
		{
			tg_selected_BCs_plot_gene_c[[i]]<-""
		}
		if(length(tg_samples_c)>0)
		{
			tg_selected_BCs_plot_sample_c[[i]]<-tg_samples_c
		}
		else
		{
			tg_selected_BCs_plot_sample_c[[i]]<-""
		}
		tg_genes_all<-c(tg_genes_all,tg_genes_c)
		tg_samples_all<-c(tg_samples_all,tg_samples_c)
		tg_col_c3<-c(tg_col_c3,rep(hmcols2[i],length(tg_samples_c)))
		tg_row_c3<-c(tg_row_c3,rep(hmcols2[i],length(tg_genes_c)))
	}
	#library(gplots)
	colors = c(-100:100)/100
	my_palette <- grDevices::colorRampPalette(c("red","white", "blue"))(n = 200)
	colors2 = c(-100:100)/100
	my_palette2 <- grDevices::colorRampPalette(c("white","white", "midnightblue"))(n = 200)
	tg_bbb0<-qubic_input_all[[ii]][tg_genes_all,tg_samples_all]
	h<-heatmap.2(tg_bbb0,Rowv=F,Colv =F,scale="none",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c3),RowSideColors=as.character(tg_row_c3),
	col=my_palette,breaks=colors,density.info="none",dendrogram="none",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_bbb1<-tg_selected_data_list2[[ii]][tg_genes_all,tg_samples_all]
	h<-heatmap.2(tg_bbb1,Rowv=F,Colv =F,scale="row",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c3),RowSideColors=as.character(tg_row_c3),
	col=my_palette,breaks=colors,density.info="none",dendrogram="none",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_r_genes<-setdiff(rownames(qubic_input_all[[ii]]),tg_genes_all)
	tg_r_samples<-setdiff(colnames(qubic_input_all[[ii]]),tg_samples_all)
	cor_r_genes<-c()
	cor_r_samples<-c()
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[i]]
		tg_samples<-tg_selected_BCs[[ii]][[2]][[i]]
		cor_r_genes<-cbind(cor_r_genes,apply(cor(t(tg_selected_data_list2[[ii]][tg_r_genes,tg_samples_all]),t(tg_selected_data_list2[[ii]][tg_genes,tg_samples_all])),1,mean))
		cor_r_samples<-cbind(cor_r_samples,apply(cor(tg_selected_data_list2[[ii]][tg_genes_all,tg_r_samples],tg_selected_data_list2[[ii]][tg_genes_all,tg_samples]),1,mean))
	}
	rownames(cor_r_samples)<-tg_r_samples
	rownames(cor_r_genes)<-tg_r_genes
	cor_r_samples_indi<-c()
	for(i in 1:nrow(cor_r_samples))
	{
		cor_r_samples_indi<-c(cor_r_samples_indi,which(cor_r_samples[i,]==max(cor_r_samples[i,]))[1])
	}
	cor_r_genes_indi<-c()
	for(i in 1:nrow(cor_r_genes))
	{
		cor_r_genes_indi<-c(cor_r_genes_indi,which(cor_r_genes[i,]==max(cor_r_genes[i,]))[1])
	}
	tg_selected_BCs_other_sample_c<-list()
	for(i in 1:length(aaa1))
	{
		if(length(names(cor_r_samples[which(cor_r_samples_indi==i),i])[order(-cor_r_samples[which(cor_r_samples_indi==i),i])])>0)
		{
			tg_selected_BCs_other_sample_c[[i]]<-names(cor_r_samples[which(cor_r_samples_indi==i),i])[order(-cor_r_samples[which(cor_r_samples_indi==i),i])]
		}
		else
		{
			tg_selected_BCs_other_sample_c[[i]]<-""
		}
	}
	tg_selected_BCs_other_gene_c<-list()
	for(i in 1:length(aaa1))
	{
		if(length(names(cor_r_genes[which(cor_r_genes_indi==i),i])[order(-cor_r_genes[which(cor_r_genes_indi==i),i])])>0)
		{
			tg_selected_BCs_other_gene_c[[i]]<-names(cor_r_genes[which(cor_r_genes_indi==i),i])[order(-cor_r_genes[which(cor_r_genes_indi==i),i])]
		}
		else
		{
			tg_selected_BCs_other_gene_c[[i]]<-""
		}
	}
	tg_col_c4<-c()
	tg_row_c4<-c()
	tg_samples_all_total<-c()
	tg_genes_all_total<-c()
	for(i in 1:length(aaa1))
	{
		tg_samples<-tg_selected_BCs[[ii]][[2]][[aaa1[i]]]
		tg_samples_c<-setdiff(tg_samples,tg_samples_all_total)
		tg_samples_all_total<-c(tg_samples_all_total,tg_samples_c)
		tg_col_c4<-c(tg_col_c4,rep(hmcols2[i],length(tg_samples_c)))
		tg_samples_all_total<-c(tg_samples_all_total,names(cor_r_samples[which(cor_r_samples_indi==aaa1[i]),aaa1[i]])[order(-cor_r_samples[which(cor_r_samples_indi==aaa1[i]),aaa1[i]])])
		tg_col_c4<-c(tg_col_c4,rep("white",length(names(cor_r_samples[which(cor_r_samples_indi==aaa1[i]),aaa1[i]])[order(-cor_r_samples[which(cor_r_samples_indi==aaa1[i]),aaa1[i]])])))	
	}
	for(i in 1:length(aaa1))
	{
		tg_genes<-tg_selected_BCs[[ii]][[1]][[aaa1[i]]]
		tg_genes_c<-setdiff(tg_genes,tg_genes_all_total)
		tg_genes_all_total<-c(tg_genes_all_total,tg_genes_c)
		tg_row_c4<-c(tg_row_c4,rep(hmcols2[i],length(tg_genes_c)))
		tg_genes_all_total<-c(tg_genes_all_total,names(cor_r_genes[which(cor_r_genes_indi==aaa1[i]),aaa1[i]])[order(-cor_r_genes[which(cor_r_genes_indi==aaa1[i]),aaa1[i]])])
		tg_row_c4<-c(tg_row_c4,rep("white",length(names(cor_r_genes[which(cor_r_genes_indi==aaa1[i]),aaa1[i]])[order(-cor_r_genes[which(cor_r_genes_indi==aaa1[i]),aaa1[i]])])))
		
	}
	tg_bbb0<-qubic_input_all[[ii]][tg_genes_all_total,tg_samples_all_total]
	h<-heatmap.2(tg_bbb0,Rowv=F,Colv =F,scale="none",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c4),RowSideColors=as.character(tg_row_c4),
	col=my_palette,breaks=colors,density.info="none",dendrogram="none",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	tg_bbb1<-tg_selected_data_list2[[ii]][tg_genes_all_total,tg_samples_all_total]
	h<-heatmap.2(tg_bbb1,Rowv=F,Colv =F,scale="row",main=tg_cancer_keys[ii],ColSideColors=as.character(tg_col_c4),RowSideColors=as.character(tg_row_c4),
	col=my_palette,breaks=colors,density.info="none",dendrogram="none",trace="none",margin=c(10,10),cexRow=0.8,cexCol=0.8)
	names(tg_selected_BCs_plot_gene_c)<-aaa1
	names(tg_selected_BCs_plot_sample_c)<-aaa1
	names(tg_selected_BCs_other_gene_c)<-1:length(aaa1)
	names(tg_selected_BCs_other_sample_c)<-1:length(aaa1)
	tg_selected_BCs_plot_c<-list(tg_selected_BCs_plot_gene_c,tg_selected_BCs_plot_sample_c)
	tg_selected_BCs_other_sample_genes_c<-list(tg_selected_BCs_other_gene_c,tg_selected_BCs_other_sample_c)
	tg_selected_BCs_plot[[ii]]<-tg_selected_BCs_plot_c
	tg_selected_BCs_other_sample_genes[[ii]]<-tg_selected_BCs_other_sample_genes_c
}
names(tg_selected_BCs_plot)<-tg_cancer_keys
names(tg_selected_BCs_other_sample_genes)<-tg_cancer_keys
return(list(tg_selected_BCs,tg_selected_BCs_plot,tg_selected_BCs_other_sample_genes))
}

Compute_Jaccard_mm<-function(M0,M1)
{
	ccc<-cbind(apply(M0,1,mean),apply(M1,1,mean))
	ddd<-sum((apply(ccc,1,min)+0.01))/sum((apply(ccc,1,max)+0.01))
	return(ddd)
}

coexpression_marker_selection2<-function(cor_C0,marker_stats1,cor_Cut_list=seq(50,90,by=2)/100)
{
	ff_list_all<-list()
	tg_markers_cor_list_all<-list()
	tg_stat_all<-c()
	for(ii in 1:ncol(marker_stats1))
	{
		print(colnames(marker_stats1)[ii])
		tg_stat_c<-rep(0,length(cor_Cut_list))
		tg_cell<-colnames(marker_stats1)[ii]
		tg_markers<-intersect(rownames(cor_C0),names(which((apply(marker_stats1!=0,1,sum)<=4)&(marker_stats1[,tg_cell]==1))))
		tg_B_markers<-tg_markers
		ff_list_c<-list()
		tg_markers_cor_list_c<-list()
		for(jj in 1:length(cor_Cut_list))
		{
			cor_Cut<-cor_Cut_list[jj]
			cn<-c()
			tg_B_markers_cor_lists<-list()
			N<-0
			for(i in 1:length(tg_B_markers))
			{
				ccc<-names(which(cor_C0[tg_B_markers[i],]>cor_Cut))
				if(length(ccc)>4)
				{
					N<-N+1
					cn<-c(cn,tg_B_markers[i])
					tg_B_markers_cor_lists[[N]]<-ccc
				}
			}
			names(tg_B_markers_cor_lists)<-cn
			ff_all<-c()
			tg_lists<-tg_B_markers_cor_lists
			if(length(tg_lists)>0)
			{
				for(i in 1:length(tg_lists))
				{
					ccc<-extract_data_symbol(marker_stats1,tg_lists[[i]])
					aaa<-unique(rownames(ccc))
					ddd<-c()
					for(j in 1:length(aaa))
					{
						if(sum(rownames(ccc)==aaa[j])==1)
						{
							bbb<-ccc[aaa[j],]
							bbb[which(bbb==0)]<-100
							bbb<-1/bbb*(bbb<=3)
							ddd<-rbind(ddd,bbb)
						}
						else
						{
							bbb<-ccc[which(rownames(ccc)==aaa[j]),]
							bbb[which(bbb==0)]<-100
							bbb<-1/bbb*(bbb<=3)
							ddd<-rbind(ddd,apply(bbb,2,mean))
						}		
					}
					rownames(ddd)<-aaa
					fff<-apply(ddd,2,mean)
					ff_all<-rbind(ff_all,fff)
				}
				if(length(ff_all)>0)
				{
					rownames(ff_all)<-names(tg_lists)
				}
			}
			tg_stat_c[jj]<-length(tg_lists)
			if(length(ff_all)==0)
			{
				ff_all<-0
				tg_B_markers_cor_lists<-0
			}
			ff_list_c[[jj]]<-ff_all
			tg_markers_cor_list_c[[jj]]<-tg_B_markers_cor_lists
		}
		names(tg_stat_c)<-cor_Cut_list
		names(ff_list_c)<-cor_Cut_list
		tg_stat_all<-cbind(tg_stat_all,tg_stat_c)
		ff_list_all[[ii]]<-ff_list_c
		tg_markers_cor_list_all[[ii]]<-tg_markers_cor_list_c
	}
	names(tg_markers_cor_list_all)<-colnames(marker_stats1)
	names(ff_list_all)<-colnames(marker_stats1)
	colnames(tg_stat_all)<-colnames(marker_stats1)
	return(list(tg_stat_all,ff_list_all,tg_markers_cor_list_all))
}


plot_cor_gene_lists<-function(data0,tg_list)
{	
	N<-0
	for(i in 1:length(tg_list))
	{
		if(sum(tg_list[[i]]=="")==0)
		{
			N<-N+1
		}
	}
	hmcols2 <- grDevices::colorRampPalette(rev(c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(N) 
	tg_genes_c<-c()
	cols<-c()
	N<-0
	for(i in 1:length(tg_list))
	{
		if(sum(tg_list[[i]]=="")==0)
		{
			N<-N+1
			tg_genes_c<-c(tg_genes_c,tg_list[[i]])
			cols<-c(cols,rep(hmcols2[N],length(tg_list[[i]])))
		}
	}
	heatmap.2(cor(t(data0[tg_genes_c,])),
	Rowv=F,Colv =F,scale="none",main=paste(tg_cancer_keys[j],names(IM_stable_genes)[i]),
	ColSideColors=as.character(cols),RowSideColors=as.character(cols),
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	trace="none",margin=c(10,10),cexRow=0.5,cexCol=1)
}

####################
analyze_QUBIC_genes_stats<-function(BC_im_depletion_gene_stats)
{
	aaa<-list()
	for(i in 1:length(BC_im_depletion_gene_stats))
	{
		aaa_c<-list()
		for(j in 1:length(BC_im_depletion_gene_stats[[i]]))
		{
			cc<-c()
			dd<-apply(BC_im_depletion_gene_stats[[i]][[j]],2,sum)[c(1,3)]
			for(k in 1:nrow(BC_im_depletion_gene_stats[[i]][[j]]))
			{	
				ff<-rbind(BC_im_depletion_gene_stats[[i]][[j]][k,c(1,3)],dd)
				pp<--1
					if((ff[1,1]/ff[2,1])>(ff[1,2]/ff[2,2]))
					{
						pp<-fisher.test(ff)$p.value
					}
				cc<-c(cc,pp)
			}
			names(cc)<-rownames(BC_im_depletion_gene_stats[[i]][[j]])
			aaa_c[[j]]<-cc
		}
		aaa[[i]]<-aaa_c
	}
	names(aaa)<-names(BC_im_depletion_gene_stats)
	return(aaa)
}


compute_Jaccard<-function(data_list)
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
				aaa[i,j]<-aaa[j,i]<-length(intersect(data_list[[i]],data_list[[j]]))/length(union(data_list[[i]],data_list[[j]]))
			}
		}
	}
	diag(aaa)<-1
	return(aaa)
}

compute_Intersect<-function(data_list)
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
				aaa[i,j]<-aaa[j,i]<-length(intersect(data_list[[i]],data_list[[j]]))
			}
		}
	}
	diag(aaa)<-0
	return(aaa)
}

norm_transform<-function(data0)
{
	data1<-data0
	for(i in 1:nrow(data0))
	{
		data1[i,]<-pnorm(data0[i,],mean(data0[i,]),sd(data0[i,]))
	}
	return(data1)
}

test_cell_list<-function(BC_gene_list,tg_cell_list,all_genes)
{
	cc_all<-list()
	for(i in 1:length(BC_gene_list))
	{
		cc<-c()
		bbb<-BC_gene_list[[i]]
		cc<-c()
		for(j in 1:length(tg_cell_list))
		{
			cc<-rbind(cc,c(length(intersect(bbb,tg_cell_list[[j]])),length(bbb),length(intersect(all_genes,tg_cell_list[[j]]))))
			#print(setdiff(bbb,all_genes))
		}
		rownames(cc)<-names(tg_cell_list)
		cc_all[[i]]<-cc
	}
	return(cc_all)
}

#########################
coexpression_marker_selection<-function(tg_markers,cor_Cut=0.8,cor_C0)
{
	tg_B_markers<-tg_markers
	cn<-c()
tg_B_markers_cor_lists<-list()
N<-0
for(i in 1:length(tg_B_markers))
{
	ccc<-names(which(cor_C0[tg_B_markers[[i]],]>cor_Cut))
	if(length(ccc)>4)
	{
		N<-N+1
		cn<-c(cn,tg_B_markers[[i]])
		tg_B_markers_cor_lists[[N]]<-ccc
	}
}
names(tg_B_markers_cor_lists)<-cn

ff_all<-c()
tg_lists<-tg_B_markers_cor_lists
for(i in 1:length(tg_lists))
{
	ccc<-extract_data_symbol(marker_stats1,tg_lists[[i]])
	aaa<-unique(rownames(ccc))
	ddd<-c()
	for(j in 1:length(aaa))
	{
		if(sum(rownames(ccc)==aaa[j])==1)
		{
			bbb<-ccc[aaa[j],]
			bbb[which(bbb==0)]<-100
			bbb<-1/bbb*(bbb<=3)
			ddd<-rbind(ddd,bbb)
		}
		else
		{
			bbb<-ccc[which(rownames(ccc)==aaa[j]),]
			bbb[which(bbb==0)]<-100
			bbb<-1/bbb*(bbb<=3)
			ddd<-rbind(ddd,apply(bbb,2,mean))
		}
	}
	rownames(ddd)<-aaa
	fff<-apply(ddd,2,mean)
	ff_all<-rbind(ff_all,fff)
}
rownames(ff_all)<-names(tg_lists)
if(nrow(ff_all)==1)
{
	ff_all<-rbind(ff_all,ff_all)	
}
if(nrow(ff_all)>1)
{
#library(gplots)
colors = c(0:100)/100
my_palette <- grDevices::colorRampPalette(c("white", "midnightblue"))(n =100)

heatmap.2(ff_all,Rowv=F,Colv =F,scale="none",main="",
	col=my_palette,breaks=colors,density.info="none",dendrogram="both",
	trace="none",margin=c(10,10),cexRow=0.6,cexCol=1)
}
	return(list(ff_all,tg_B_markers_cor_lists))
}