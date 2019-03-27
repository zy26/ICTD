
########################## function directory
#do_regression_unique
#extract_marker
#extract_marker_ICTD
#filter_gene_name
#marker_zero_length
#plot_cor
#RMSE_one_vector
#RMSE_two_vector
#R2_two_vector
#R2_two_vector_withIntercept
#RMSE_one_mat
#RMSE_two_mat
#R2_two_mat
#R2_two_mat_withIntercept
#rm_zero_row
#getFractions.Abbas
#TIMER_TCGA
#TIMER_sc
#unique_signature
#validate_marker_coverage_by_proportion
#validate_marker_std_output_parameter_4.6_4.7	#NMFinput may report bug (matrix of list)
#############################

 do_regression_unique <- function(proportion=Prop_EPIC[,1:7], data=bulk1, indicate=ss_signature){
 	print("unique regression")
 	if(nrow(data) != nrow(indicate)) stop("number of gene in indicate NOT match!")
 	if(ncol(proportion) != ncol(indicate)) stop("number of sample in indicate NOT match!")
 	P_2nd <- matrix(0, nrow(indicate), ncol(indicate))
 	rownames(P_2nd) <- rownames(indicate)
	n_gene <- nrow(data)
	n_sample <- ncol(data)	
	for(i in 1:n_gene){
		#print(i)
		cell <- which(indicate[i, ] == 1)
		#if(length(cell) > 1) stop("return 2 cell type!")
	 	pp <- proportion[, cell]
		y <- data[i, ]
		ddd <- cbind(pp, y)
		ddd <- as.data.frame(ddd)
		coeff <- lm(y~pp + 0, ddd)	#NOTE : the order is matter!
		P_2nd[i, cell] <- coeff[[1]]
	}	
	return(P_2nd)
 }


extract_marker <- function(aaa, method, add){
	# mmmm is a list to store cell marker, the name of each element is cell type,
	# 
	cutoff <- 0.05
	mmm <- matrix(0, nrow(aaa), ncol(aaa))
	rownames(mmm) <- rownames(aaa)
	colnames(mmm) <- colnames(aaa)
	for(i in 1:nrow(aaa)){
		vv <- matrix(NA, 1, ncol(aaa))
		vv <- aaa[i, ]

		gene_sd <- sd(vv) * sqrt( (length(vv)-1) / (length(vv))  )
		gene_mean <- mean(vv)

		#z-score calculation
		z <- (vv - gene_mean) / gene_sd

		p_yellow <- pnorm(z)
		p_blue <- 1 - p_yellow
		#which(p_blue < cutoff)
		mmm[i, which(p_blue < cutoff)] <- 1
	}

	#find zero rows
	#mmm <- extract_marker(aaa)
	#row_sub <- apply(mmm, 1, function(row) all( row == 0)) #return logistic value(T/F)
	#Zero <- which(row_sub == T)
	#length(Zero) 
	#length(which(mmm[, 1] != 0))
	#length(which(mmm[, 2] != 0) )
	#length(which(mmm[, 3] != 0) )
	#length(which(mmm[, 4] != 0) )
	#length(which(mmm[, 5] != 0) )
	#length(which(mmm[, 6] != 0) )

	if(method == "TIMER"){
		B_mark <- rownames(mmm)[which(mmm[, 1] != 0) ]
		CD4_mark <- rownames(mmm)[which(mmm[, 2] != 0) ]
		CD8_mark <-rownames(mmm)[which(mmm[, 3] != 0) ]
		Netro_mark <- rownames(mmm)[which(mmm[, 4] != 0) ]
		Macro_mark <- rownames(mmm)[which(mmm[, 5] != 0) ]
		Dc_mark <- rownames(mmm)[which(mmm[, 6] != 0) ]
	}
	#epic
	if(method == "EPIC"){
		B_mark <- rownames(mmm)[which(mmm[, 1] != 0) ]
		CAFs_mark <- rownames(mmm)[which(mmm[, 2] != 0) ]
		CD4_mark <-rownames(mmm)[which(mmm[, 3] != 0) ]
		CD8_mark <- rownames(mmm)[which(mmm[, 4] != 0) ]
		endoth_mark <- rownames(mmm)[which(mmm[, 5] != 0) ]
		Macro_mark <- rownames(mmm)[which(mmm[, 6] != 0) ]
		NK_mark <- rownames(mmm)[which(mmm[, 7] != 0) ]

		if(add==T) CD8_mark <- c(CD8_mark, "TSPYL1","UBASH3A")

	}
	#ciber
	if(method == "CIBER"){
		specific_marker <- list()
		for(i in 1:ncol(mmm)){
			specific_marker[[i]] <- rownames(mmm)[which(mmm[, i] != 0) ]
			names(specific_marker)[i] <- colnames(mmm)[i]
		}
	}


	#return
	if(method == "TIMER"){
		return(list(B_mark = B_mark, CD4_mark = CD4_mark, CD8_mark = CD8_mark, Netro_mark = Netro_mark, Macro_mark =Macro_mark, DC_mark = Dc_mark))
	}
	if(method == "EPIC"){
		return(list(B_mark = B_mark, CAFs_mark = CAFs_mark, CD4_mark = CD4_mark, CD8_mark = CD8_mark, endoth_mark =endoth_mark, Macro_mark = Macro_mark, NK_mark=NK_mark))
	}
	if(method == "CIBER"){
		return(specific_marker)
	}

}


extract_marker_ICTD <- function(C){
	specific_marker <- list()
	for(i in 1:ncol(C)){
		specific_marker[[i]] <- rownames(C)[which(C[, i] != 0)]
		names(specific_marker)[i] <- colnames(C)[i]
	}

	return(specific_marker)
}


filter_gene_name <- function(data_t){
# 1. extra special row
loca <- c()
for(i in 1:nrow(data_t)){
	tmp <- unlist(strsplit(rownames(data_t)[i], '|', fixed = T))
	if(length(tmp) > 1) loca <- c(loca, i)
}
data_sub <- data_t[loca,]

# 2. change name
gene <- c()
for(i in 1:nrow(data_sub)){
	tmp <- unlist(strsplit(rownames(data_sub)[i], '|', fixed = T))
	gene <- c(gene, tmp)
}
gene_2col <- matrix(gene, 2, length(gene)/2)
gene_name <- gene_2col[2,]
rownames(data_sub) <- gene_name

data_ttt <- data_sub

return(data_ttt)

}

marker_zero_length <- function(aaaa){
	zero_length = 0
	for(i in 1:length(aaaa)){
		if(length(aaaa[[i]]) == 0) zero_length=zero_length+1
	}

	return(zero_length)
}

RMSE_one_vector <- function(a){
	rmse <- sqrt(sum(a^2)/length(a))
	return(rmse)
}

RMSE_two_vector <- function(a, b){
	diff <- a - b
	rmse <- sqrt(sum(diff^2)/length(diff))
	return(rmse)
}

R2_two_vector <- function(predict, actual){
	#R2 <- 1 - sum( (actual-predict )^2 ) / sum( (actual-mean(actual) )^2  ) 
	R2 <- 1 - sum( (actual-predict )^2 ) / sum( actual^2 ) 
	return(R2)
}

R2_two_vector_withIntercept <- function(predict, actual){
	R2 <- 1 - sum( (actual-predict )^2 ) / sum( (actual-mean(actual) )^2  ) 
	return(R2)
}

RMSE_one_mat <- function(aaa){
	
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		vc[i] <- RMSE_one_vector(tmp_a)

	}
	rownames(vc) <- rownames(aaa)
	
	return(vc)
}

RMSE_two_mat <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- RMSE_two_vector(tmp_a, tmp_b)

	}
	rownames(vc) <- rownames(aaa)
	
	return(vc)
}

R2_two_mat <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- R2_two_vector(tmp_a, tmp_b)
	}
	rownames(vc) <- rownames(aaa)
	return(vc)
}

R2_two_mat_withIntercept <- function(aaa, bbb){
	if(nrow(aaa) != nrow(bbb)) stop("size of aaa and bbb different!")
	n_gene <- nrow(aaa)
	vc <- matrix(NA, n_gene, 1)
	for(i in 1:n_gene){
		tmp_a <- aaa[i, ]
		tmp_b <- bbb[i, ]
		vc[i] <- R2_two_vector_withIntercept(tmp_a, tmp_b)
	}
	rownames(vc) <- rownames(aaa)
	return(vc)
}


rm_zero_row <- function(bulk){

	row_sub <- apply(bulk, 1, function(row) all( row == 0)) #return logistic value(T/F)
	Zero <- which(row_sub == TRUE)
	print(length(Zero))
	if(length(Zero) == 0) bulk <- bulk
	if(length(Zero) > 0)  bulk <- bulk[-Zero,]

	return(bulk)
}

getFractions.Abbas <- function(XX,YY,w=NA){
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0,ncol(XX)))
      tmp.XX=tmp.XX[,-ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0,ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX,YY,intercept=F) else tmp=lsfit(tmp.XX,YY,w,intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0,ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

TIMER_TCGA <- function(cc = cancer_str){
	cc = tolower(cc)
	if(cc=='skcm')cc.type='06A' else cc.type='01A'

	##----- setup parameters and establish the output file -----##
	signature.genes=c('CD19','TRAT1','CD8B','CCR3','CD163','CCL17')
	names(signature.genes)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

	##----- load and process gene expression data -----##
	file_str <- paste("C:/Users/wnchang/Documents/F/PhD_Research/TCGA_data/TCGA-", toupper(cancer_str), "_FPKM_T.RData", sep="")
	data_t <- get(load(file_str))
	data_ttt <- filter_gene_name(data_t)
	dd <- data_ttt
	#dd <- get(load(file_str) )
	dd=as.matrix(dd)
	mode(dd)='numeric'	
	#if(!cc %in% c('gbm','ov','esca','stad'))dd=dd*1e6   ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
	tmp=strsplit(rownames(dd),'\\|')
	tmp=sapply(tmp,function(x)x[[1]])
	tmp.vv=which(nchar(tmp)>1)
	rownames(dd)=tmp
	dd=dd[tmp.vv,]

	##----- load immune marker genes from Abbas et al., 2005 -----##
	tmp=read.csv("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/IRIS-marker-gene.txt", ,header=T,sep='\t',stringsAsFactors=F) # * add / before data
	marker.list=tmp[,1]
	names(marker.list)=tmp[,7]
	names(marker.list)=gsub(' ','_',tmp[,7])
	names(marker.list)=gsub('Dendritic_Cell','DC',names(marker.list))
	names(marker.list)=gsub('Neutrophil','Neutrophils',names(marker.list))
	gsub('Cell','cell',names(marker.list))->names(marker.list)

	##----- load reference data of sorted immune cells -----##
	load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/HPCTimmune.Rdata") # * add / before data

	##----- load and process tumor purity data -----##
	AGP=read.table(paste('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/AGP/AGP-',cc,'.txt',sep=''),sep='\t',header=T) # * add / before datas
	AGP=AGP[which(AGP[,'PoP']>0.01),]
	tmp=strsplit(rownames(HPCT.immune),';')
	AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
	names(AffyIDtoGenes)=sapply(tmp,function(x)x[[2]])
	marker.list.genes=AffyIDtoGenes[marker.list]

	##----- function to edit TCGA ID, with the option of keeping the first num.res fields -----##
	getID <- function(sID,num.res=3){
 		mm=c()
  		for(id in sID){
   			tmp=unlist(strsplit(id,'-'))
    		if(length(tmp)==1){
     	 		tmp=unlist(strsplit(id,'\\.'))
    		}
    		ll='TCGA'
    		for(j in 2:num.res){
     			ll=paste(ll,tmp[j],sep='-')
   		 	}
    		mm=c(mm,ll)
 		}
  		return(mm)
	}
	rownames(AGP)=getID(AGP[,1],4)
	colnames(dd)=getID(colnames(dd),4)


	##----- Select single reference samples of pre-selected immune cell types -----##
	B_cell=362:385
	T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
	T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
	NK=328:331
	Neutrophil=344:361
	Macrophage=66:80
	DC=151:238
	curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

	curated.cell.types=colnames(curated.ref)
	names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))

	load('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/curated.ref.genes.Rdata')

	##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
	RemoveBatchEffect <- function(){
  		library(sva)
  		tmp.dd=as.matrix(dd)
  		tmp=sapply(strsplit(rownames(dd),'\\|'),function(x)x[[1]])
  		rownames(tmp.dd)=tmp
  		tmp.dd=tmp.dd[which(nchar(tmp)>1),]
  		tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  		N1=ncol(tmp.dd)
  		tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,])
  		tmp.dd=as.matrix(tmp.dd)
  		mode(tmp.dd)='numeric'
  		N2=ncol(curated.ref.genes)
  		tmp.batch=c(rep(1,N1),rep(2,N2))
  		tmp.dd0=ComBat(tmp.dd,tmp.batch,c())
  		dd.br=tmp.dd0[,1:N1]
  		curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)]
  		tmp0=c()
  		for(kk in unique(names(curated.cell.types))){
   		 tmp.vv=which(names(curated.cell.types)==kk)
   		 tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  		}
  		curated.ref.genes.agg.br=tmp0
  		colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  		#rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  		return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
	}

	tmp=RemoveBatchEffect()
	dd.br=tmp$dd
	curated.ref.genes.br=tmp$rr
	curated.ref.genes.agg.br=tmp$rrg


	##----- function to calculate the residuals from regression -----##
	fn <- function(beta0,XX,Y)return(log(sum(abs(Y-XX%*%beta0))))

	##----- function to select genes with expression values negatively correlated with tumor purity -----##
	getPurityGenes <- function(dd,AGP,thr.p=0.05,thr.c=0,mode='env'){
 	 tmp.ss=intersect(colnames(dd),rownames(AGP))
  	if(length(tmp.ss)==0){
    	colnames(dd)=getID(colnames(dd))
   		tmp.ss=intersect(colnames(dd),rownames(AGP))
  	}
  	tmp.dd=dd[,tmp.ss]	#original
 	#tmp.dd <- dd    #w
  	tmp=lapply(rownames(tmp.dd),function(x)cor.test(tmp.dd[x,],as.numeric(AGP[colnames(tmp.dd),2]),method='s'))
  	tmp.pp=sapply(tmp,function(x)x$p.value)
  	tmp.cor=sapply(tmp,function(x)x$estimate)
  	names(tmp.pp)=names(tmp.cor)=rownames(dd)
  	if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
 	if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
 	return(vv)
	}

	##----- selection genes negatively correlated with purity and overlap with immune marker genes -----##
	vv.t=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= -0.2)
	vv.t=intersect(vv.t,rownames(curated.ref.genes.agg.br))
	vv=intersect(vv.t,marker.list.genes)

	##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
	RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
 	 ## removes upper thr.q quantile for every reference feature
 	 remove.vv=c()
 	 for(i in 1:ncol(ref.dd)){
  	  tmp=quantile(ref.dd[vv,i],thr.q)[1]
  	  tmp.vv=which(ref.dd[vv,i]>tmp)
 	  remove.vv=c(remove.vv,tmp.vv)
	 }
 	 remove.vv=unique(remove.vv)
 	 return(vv[-remove.vv])
	}

	##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
	tmp.diff=sum(sum(abs(cor(curated.ref.genes.agg.br[vv,],method='p')-cor(curated.ref.genes.agg.br[vv,],method='s'))))

	if(tmp.diff>= -10000){
 		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4])
 		 vv=vv0
	}

	cat("Number of genes inversely correlated with purity is ",length(vv.t),'\n\n',sep='',file='output-statistics.txt',append=T)
	cat("Number of immune genes inversely correlated with purity is ",length(vv),'\n\n',sep='',file='output-statistics.txt',append=T)

	##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
	tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
	n.immune=length(intersect(marker.list.genes,tmp.ss0))
	cat("Test if immune genes are enriched for inverse correlation with purity: \n\n",file='output-statistics.txt',append=T)
	sink(file='output-statistics.txt',append=T);print(fisher.test(matrix(c(length(vv),length(vv.t)-length(vv),n.immune,length(tmp.ss0)-n.immune),2,2)));sink()


	##----- function to process deconvolution method in batch -----##
	BatchFractions <- function(XX,YYd){
  		Fmat=c()
  		for(i in 1:ncol(YYd)){
  			YY=YYd[,i]
   			tmp.F=getFractions.Abbas(XX,YY)
    		#tmp.F=getFractions.Optim(XX,YY)
    		Fmat=rbind(Fmat,tmp.F)
  		}
  		rownames(Fmat)=colnames(YYd)
  		colnames(Fmat)=colnames(XX)
  		return(Fmat)
	}

	##----- perform batch deconvolution -----##
	XX=curated.ref.genes.agg.br[vv,c(-4)]	#chang note: delete NK cell during deconvolution
	YYd=dd.br[vv,]
	Fmat=BatchFractions(XX,YYd)


	##----- CD4 and CD8 T cells are likely to be similar, resulting in colinearity. Codes below are procedures to remove outlier genes that may result in colinearity until the two covariates are linearly separable. -----##
	if(cor(Fmat[,2],Fmat[,3])<= -0.2){
  		if(tmp.diff>=1){
    		tmp.cor=c()
    		thr.qlist=c(0.99)
    		for(tq in thr.qlist){
     		 vv=intersect(vv.t,marker.list.genes)
    		 vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4],tq)
    		 vv=vv0
     		 XX=curated.ref.genes.agg.br[vv,]
     		 YYd=dd.br[vv,]
     		 tmp.Fmat=BatchFractions(XX,YYd)
     		 tmp.cor=c(tmp.cor,cor(tmp.Fmat[,2],tmp.Fmat[,3],method='s'))
   			}
   			tmp.vv=which.max(tmp.cor)
    		tq=thr.qlist[tmp.vv]
    		vv=intersect(vv.t,marker.list.genes)
    		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,c(-4)],tq)
    		vv=vv0
    		XX=curated.ref.genes.agg.br[vv,c(-4)]
    		YYd=dd.br[vv,]
    		Fmat=BatchFractions(XX,YYd)
    		Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
    		rownames(Fmat0.p)=getID(rownames(Fmat0.p))
  		}
	}

	while(cor(Fmat[,2],Fmat[,3])<=-0.3){
 	 if(length(vv)<=50)break
	 vv=vv[-as.numeric(names(table(apply(dd[vv,],2,which.max))))]
  	 XX=curated.ref.genes.agg.br[vv,c(-4)]
     #XX=XX[,!colnames(XX) %in% tmp.remove]
     YYd=dd.br[vv,]
     Fmat=BatchFractions(XX,YYd)
	}

	Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
	rownames(Fmat0.p)=getID(rownames(Fmat0.p))

	TIMER_fraction <- Fmat
	X <- YYd
	S <- XX #negative value
	return(list(Timer_fraction = TIMER_fraction, X = X, S = S))

}

TIMER_sc <- function(cc = cancer_str, data.matrix = data.matrix){
	cc = tolower(cc)
	cancer_lib <- c("BLCA", "BRCA", "CESC", "COAD", "DLBC","ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", 
				"LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "UCS", "UVM")
	cancer_lib_lll <- tolower(cancer_lib)
	if(!cc%in%cancer_lib_lll)
		stop("please input the cancer type which exist in TCGA library. For example BLCA, COAD, SKCM...")
	if(cc=='skcm')cc.type='06A' else cc.type='01A'

	##----- setup parameters and establish the output file -----##
	signature.genes=c('CD19','TRAT1','CD8B','CCR3','CD163','CCL17')
	names(signature.genes)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

	##----- load and process gene expression data -----##
	#check sample name
	if(length(colnames(data.matrix)) == 0) {
	warning("input data do NOT have colnames")
	colnames(data.matrix) <- paste( "Setsample", 1:ncol(data.matrix), sep="")    
	}
	data.matrix <- rm_zero_row(data.matrix)
	#data_ttt <- filter_gene_name(data_t)
	data_ttt <- data.matrix
	dd <- data_ttt
	#dd <- get(load(file_str) )
	dd=as.matrix(dd)
	mode(dd)='numeric'	
	#if(!cc %in% c('gbm','ov','esca','stad'))dd=dd*1e6   ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
	tmp=strsplit(rownames(dd),'\\|')
	tmp=sapply(tmp,function(x)x[[1]])
	tmp.vv=which(nchar(tmp)>1)
	rownames(dd)=tmp
	dd=dd[tmp.vv,]

	##----- load immune marker genes from Abbas et al., 2005 -----##
	tmp=read.csv("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/IRIS-marker-gene.txt", ,header=T,sep='\t',stringsAsFactors=F) # * add / before data
	marker.list=tmp[,1]
	names(marker.list)=tmp[,7]
	names(marker.list)=gsub(' ','_',tmp[,7])
	names(marker.list)=gsub('Dendritic_Cell','DC',names(marker.list))
	names(marker.list)=gsub('Neutrophil','Neutrophils',names(marker.list))
	gsub('Cell','cell',names(marker.list))->names(marker.list)

	##----- load reference data of sorted immune cells -----##
	load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/HPCTimmune.Rdata") # * add / before data

	##----- load and process tumor purity data -----##
	AGP=read.table(paste('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/AGP/AGP-',cc,'.txt',sep=''),sep='\t',header=T) # * add / before datas
	AGP=AGP[which(AGP[,'PoP']>0.01),]
	tmp=strsplit(rownames(HPCT.immune),';')
	AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
	names(AffyIDtoGenes)=sapply(tmp,function(x)x[[2]])
	marker.list.genes=AffyIDtoGenes[marker.list]

	##----- function to edit TCGA ID, with the option of keeping the first num.res fields -----##
	getID <- function(sID,num.res=3){
 		mm=c()
  		for(id in sID){
   			tmp=unlist(strsplit(id,'-'))
    		if(length(tmp)==1){
     	 		tmp=unlist(strsplit(id,'\\.'))
    		}
    		ll='TCGA'
    		for(j in 2:num.res){
     			ll=paste(ll,tmp[j],sep='-')
   		 	}
    		mm=c(mm,ll)
 		}
  		return(mm)
	}
	rownames(AGP)=getID(AGP[,1],4)
	#colnames(dd)=getID(colnames(dd),4)


	##----- Select single reference samples of pre-selected immune cell types -----##
	B_cell=362:385
	T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
	T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
	NK=328:331
	Neutrophil=344:361
	Macrophage=66:80
	DC=151:238
	curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

	curated.cell.types=colnames(curated.ref)
	names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))

	load('C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/immune_datasets/curated.ref.genes.Rdata')

	##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
	RemoveBatchEffect <- function(){
  		library(sva)
  		tmp.dd=as.matrix(dd)
  		tmp=sapply(strsplit(rownames(dd),'\\|'),function(x)x[[1]])
  		rownames(tmp.dd)=tmp
  		tmp.dd=tmp.dd[which(nchar(tmp)>1),]
  		tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  		N1=ncol(tmp.dd)
  		tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,])
  		tmp.dd=as.matrix(tmp.dd)
  		mode(tmp.dd)='numeric'
  		N2=ncol(curated.ref.genes)
  		tmp.batch=c(rep(1,N1),rep(2,N2))
  		tmp.dd0=ComBat(tmp.dd,tmp.batch,c())
  		dd.br=tmp.dd0[,1:N1]
  		curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)]
  		tmp0=c()
  		for(kk in unique(names(curated.cell.types))){
   		 tmp.vv=which(names(curated.cell.types)==kk)
   		 tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  		}
  		curated.ref.genes.agg.br=tmp0
  		colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  		#rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  		return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
	}

	tmp=RemoveBatchEffect()
	dd.br=tmp$dd
	curated.ref.genes.br=tmp$rr
	curated.ref.genes.agg.br=tmp$rrg


	##----- function to calculate the residuals from regression -----##
	fn <- function(beta0,XX,Y)return(log(sum(abs(Y-XX%*%beta0))))

	##----- function to select genes with expression values negatively correlated with tumor purity -----##
	getPurityGenes <- function(dd,AGP,thr.p=0.05,thr.c=0,mode='env'){
 #	 tmp.ss=intersect(colnames(dd),rownames(AGP))
 # 	if(length(tmp.ss)==0){
 #   	colnames(dd)=getID(colnames(dd))
 #  		tmp.ss=intersect(colnames(dd),rownames(AGP))
 # 	}
 # 	tmp.dd=dd[,tmp.ss]	#original
 # 	#tmp.dd <- dd    #w
 # 	tmp=lapply(rownames(tmp.dd),function(x)cor.test(tmp.dd[x,],as.numeric(AGP[colnames(tmp.dd),2]),method='s'))
 # 	tmp.pp=sapply(tmp,function(x)x$p.value)
 # 	tmp.cor=sapply(tmp,function(x)x$estimate)
 # 	names(tmp.pp)=names(tmp.cor)=rownames(dd)
 # 	if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
# 	if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
 #	return(vv)

 		#single cell data do NOT have overlap sample with AGP, thus we assume all genes are purity genes.(chang.)
 		return(rownames(dd))
	}

	##----- selection genes negatively correlated with purity and overlap with immune marker genes -----##
	vv.t=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= -0.2)
	vv.t=intersect(vv.t,rownames(curated.ref.genes.agg.br))
	vv=intersect(vv.t,marker.list.genes)

	##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
	RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
 	 ## removes upper thr.q quantile for every reference feature
 	 remove.vv=c()
 	 for(i in 1:ncol(ref.dd)){
  	  tmp=quantile(ref.dd[vv,i],thr.q)[1]
  	  tmp.vv=which(ref.dd[vv,i]>tmp)
 	  remove.vv=c(remove.vv,tmp.vv)
	 }
 	 remove.vv=unique(remove.vv)
 	 return(vv[-remove.vv])
	}

	##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
	tmp.diff=sum(sum(abs(cor(curated.ref.genes.agg.br[vv,],method='p')-cor(curated.ref.genes.agg.br[vv,],method='s'))))

	if(tmp.diff>= -10000){
 		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4])
 		 vv=vv0
	}

	cat("Number of genes inversely correlated with purity is ",length(vv.t),'\n\n',sep='',file='output-statistics.txt',append=T)
	cat("Number of immune genes inversely correlated with purity is ",length(vv),'\n\n',sep='',file='output-statistics.txt',append=T)

	##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
	tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
	n.immune=length(intersect(marker.list.genes,tmp.ss0))
	cat("Test if immune genes are enriched for inverse correlation with purity: \n\n",file='output-statistics.txt',append=T)
	sink(file='output-statistics.txt',append=T);print(fisher.test(matrix(c(length(vv),length(vv.t)-length(vv),n.immune,length(tmp.ss0)-n.immune),2,2)));sink()


	##----- function to process deconvolution method in batch -----##
	BatchFractions <- function(XX,YYd){
  		Fmat=c()
  		for(i in 1:ncol(YYd)){
  			YY=YYd[,i]
   			tmp.F=getFractions.Abbas(XX,YY)
    		#tmp.F=getFractions.Optim(XX,YY)
    		Fmat=rbind(Fmat,tmp.F)
  		}
  		rownames(Fmat)=colnames(YYd)
  		colnames(Fmat)=colnames(XX)
  		return(Fmat)
	}

	##----- perform batch deconvolution -----##
	XX=curated.ref.genes.agg.br[vv,c(-4)]	#chang note: delete NK cell during deconvolution
	YYd=dd.br[vv,]
	Fmat=BatchFractions(XX,YYd)

if( !is.na( cor(Fmat[,2],Fmat[,3])[1] ) ){
	print("Remove co-linear!!!")
	##----- CD4 and CD8 T cells are likely to be similar, resulting in colinearity. Codes below are procedures to remove outlier genes that may result in colinearity until the two covariates are linearly separable. -----##
	if(cor(Fmat[,2],Fmat[,3])<= -0.2){
  		if(tmp.diff>=1){
    		tmp.cor=c()
    		thr.qlist=c(0.99)
    		for(tq in thr.qlist){
     		 vv=intersect(vv.t,marker.list.genes)
    		 vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4],tq)
    		 vv=vv0
     		 XX=curated.ref.genes.agg.br[vv,]
     		 YYd=dd.br[vv,]
     		 tmp.Fmat=BatchFractions(XX,YYd)
     		 tmp.cor=c(tmp.cor,cor(tmp.Fmat[,2],tmp.Fmat[,3],method='s'))
   			}
   			tmp.vv=which.max(tmp.cor)
    		tq=thr.qlist[tmp.vv]
    		vv=intersect(vv.t,marker.list.genes)
    		vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,c(-4)],tq)
    		vv=vv0
    		XX=curated.ref.genes.agg.br[vv,c(-4)]
    		YYd=dd.br[vv,]
    		Fmat=BatchFractions(XX,YYd)
    		Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
    		rownames(Fmat0.p)=getID(rownames(Fmat0.p))
  		}
	}

	count=0
	while(cor(Fmat[,2],Fmat[,3])<=-0.3){
		print("Remove co-linearality iterativelly!!~")
		count=count+1
		#print(count)
 	 if(length(vv)<=50)break
	 vv=vv[-as.numeric(names(table(apply(dd[vv,],2,which.max))))]
  	 XX=curated.ref.genes.agg.br[vv,c(-4)]
     #XX=XX[,!colnames(XX) %in% tmp.remove]
     YYd=dd.br[vv,]
     Fmat=BatchFractions(XX,YYd)
     	if(count>=30)break
	}
}
	#Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
	#rownames(Fmat0.p)=getID(rownames(Fmat0.p))

	TIMER_fraction <- Fmat
	X <- YYd
	S <- XX #negative value
	return(list(Timer_fraction = TIMER_fraction, X = X, S = S))

}

unique_signature <- function(aaa, cutoff=0.05){

	#cutoff <- 0.05
	mmm <- matrix(0, nrow(aaa), ncol(aaa))
	rownames(mmm) <- rownames(aaa)
	colnames(mmm) <- colnames(aaa)
	for(i in 1:nrow(aaa)){
		vv <- matrix(NA, 1, ncol(aaa))
		vv <- aaa[i, ]

		gene_sd <- sd(vv) * sqrt( (length(vv)-1) / (length(vv))  )
		gene_mean <- mean(vv)

		#z-score calculation
		z <- (vv - gene_mean) / gene_sd

		p_yellow <- pnorm(z)
		p_blue <- 1 - p_yellow
		#which(p_blue < cutoff)
		mmm[i, which(p_blue < cutoff)] <- 1
	}

	return(mmm)
}

#e.g.:
#overlap_gene = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],table=immune_cell_uni_table0_GS,NMFmarker=NMF_input1[[2]],cor.cutoff=0.85)
#overlap_gene4 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],table=immune_cell_uni_table0_GS,cor.cutoff=0.85)
validate_marker_coverage_by_proportion <- function(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff=0.9){
	table.gene = rownames(table)
	origiData_01 = origiData[intersect(table.gene,rownames(origiData)), ]

	ccc = cor(t(origiData_01), t(truePro))		#ccc_overlap
	ccc.whole = cor(t(origiData),t(truePro))
	ccc.order = list()							#ccc_overlap_order
	ccc.whole.order = list()

	for(i in 1:ncol(ccc)){
		ccc.order[[i]]= ccc[order(-ccc[,i])[1:100], i]
	}
	if(length(rownames(tProp)) == 0) stop("tProp do NOT has rowname(cell type name)")
	names(ccc.order) = rownames(tProp)
	for(i in 1:ncol(ccc.whole)){
		ccc.whole.order[[i]] = ccc.whole[order(-ccc.whole[,i])[1:100], i]
	}
	names(ccc.whole.order) = rownames(tProp)

	marker_true = list()
	# 0.1 version
	#for(i in 1:ncol(ccc)){	#####!!!omit 0.9 cutoff
	#	marker_true[[i]]=names(sort(ccc[,i], decreasing=TRUE)[1:100])
	#}
	#names(marker_true) = colnames(ccc)
	# 0.2 version
	for(i in 1:ncol(ccc)){	####omit top 100 limit, make sense
		cor.sort = sort(ccc[,i], decreasing=TRUE)
		marker_true[[i]] = names(cor.sort[which(cor.sort> cor.cutoff)])
	}
	names(marker_true) = colnames(ccc)

	#print true marker information
	#print("marker_true_total:")
	for(i in 1:ncol(ccc)){
		show = paste(names(marker_true)[i],"::::",length(marker_true[[i]]),sep="")
		#print(show)

	}

	#check high correlation marker do not overlap between cell type
	for(i in 1:length(marker_true)){
		for(j in (i+1):length(marker_true)){
			if(j > length(marker_true)) break
			if(length(intersect(marker_true[[i]], marker_true[[j]]))>0) {
				warning(paste(names(marker_true)[i],"--",names(marker_true)[j],"__overlap_gene:",length(intersect(marker_true[[i]], marker_true[[j]])),sep=""))
				print(intersect(marker_true[[i]], marker_true[[j]]))
			}
		}
	}

	common_mat = matrix(NA, ncol(ccc)+2, length(markList))
	for(i in 1:length(markList)){
		common_mat[ncol(ccc)+1, i] = length(markList[[i]])
	}

	common_NMF = matrix(NA, ncol(ccc)+2, ncol(NMFmarker))
	for(i in 1:ncol(NMFmarker)){
		common_NMF[ncol(ccc)+1, i] = length(which(NMFmarker[,i] == 1))
	}

	marker_miss = list()
	#print("marker_remain:")
	for(i in 1:length(marker_true)){
		geneLib = marker_true[[i]] 	#cell type in true marker list
		geneLib1 = marker_true[[i]]	#for matrix
		
		for(pp in 1:length(markList)){
			#print(names(markList)[pp])
			#print("::::::::::::::")
			geneSet = markList[[pp]]
			geneLib = setdiff(geneLib, geneSet)
			#print(length(geneLib) )
			common_mat[i,pp]=length(intersect(geneSet,geneLib1))
			
		}

		marker_miss[[i]] = geneLib
		show = paste(names(marker_true)[i],"---",length(geneLib),sep="")
		#print(show)
	}
	names(marker_miss)=names(marker_true)	#this is residual gene

	vv.tmp = apply(common_mat[1:ncol(ccc),], 2, max)
	common_mat[ncol(ccc)+2, ] = vv.tmp / common_mat[ncol(ccc)+1, ]	
	rownames(common_mat) = c(colnames(ccc),"model_size","top1_percentage")
	colnames(common_mat) = names(markList)	

	for(i in 1:length(marker_true)){
		geneLib = marker_true[[i]]

		for(k in 1:ncol(NMFmarker)){
			geneSet = names(which(NMFmarker[, k]==1) )
			common_NMF[i, k] = length(intersect(geneSet, geneLib))
		}
	}
	nn.tmp = apply(common_NMF[1:ncol(ccc), ], 2, max)
	common_NMF[ncol(ccc)+2, ] = nn.tmp / common_NMF[ncol(ccc)+1, ]
	rownames(common_NMF) = c(colnames(ccc),"model_size","top1_percentage")
	colnames(common_NMF) = colnames(NMFmarker)

	#top1 row basis (gene*) correlated with all gene in original data
	geneX = c()
	count = 0
	for(i in 1:ncol(NMFmarker)){
		gene_row= names(which(NMFmarker[,i]==1))
		gene_row_common = intersect(gene_row,rownames(origiData))
		if(length(gene_row) != 0){
			data_row = origiData[gene_row_common,]
			vv = (svd(data_row)$v)[,1]
			if(sum(vv) < 0) vv = -vv
			geneX = rbind(geneX, vv)
			count = count + 1
		}
	}
	rownames(geneX) = paste("geneX_",colnames(NMFmarker)[1:count],sep="")
	cor_top_row = cor(t(tProp), t(geneX))

	return(list(common_mat=common_mat,marker_remain=marker_miss,marker_true=marker_true,ccc_overlap=ccc,ccc_whole=ccc.whole,
		ccc_overlap_order=ccc.order,ccc_whole_order=ccc.whole.order,common_NMF=common_NMF,cor_top_row=cor_top_row))
}


validata_R4 <- function(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]]){

	geneX=c()
	for(i in 1:length(markList)){
		gene_row = markList[[i]]
		data_row = origiData[gene_row,]
		vv = (svd(data_row)$v)[, 1]
		if(sum(vv) < 0) vv = -vv
		geneX = rbind(geneX, vv)
	}
	rownames(geneX) = paste("geneX_", names(markList), "~",c(1:length(markList)), sep="")
	cor_top_row = cor(t(truePro),t(geneX))

	return(cor_top_row)
}


validate_marker_std_output <- function(data.matrix=data.matrix,tProp=tProp,R1_filter_step1_results_new=R1_filter_step1_results_new,NMF_input1=NMF_input1,immune_cell_uni_table0_GS=immune_cell_uni_table0_GS,cor.cutoff=0.9,comb=comb){

	output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)
	output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)

	o1_top100_ICTDmarker = output2[["ccc_overlap_order"]]
	o1_top100_cor_gene = output2[["ccc_whole_order"]]
	o2_R1_overlap = output1[["common_mat"]]
	o3_R4_overlap = output2[["common_mat"]]
	o3_R4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]])
	o4_R4_4.1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=comb)
	o4_R4_4.2_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=comb)
	o4_R4_4.3_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=comb)
	o5_nmf_overlap = output2[["common_NMF"]]
	o5_nmf_svd_cor = output2[["cor_top_row"]]

	return(list(o1_top100_ICTDmarker=o1_top100_ICTDmarker,o1_top100_cor_gene=o1_top100_cor_gene,o2_R1_overlap=o2_R1_overlap,
				o3_R4_overlap=o3_R4_overlap,o3_R4_svd_cor=o3_R4_svd_cor,
				o4_R4_step2_svd_cor=o4_R4_step2_svd_cor,
				o5_nmf_overlap=o5_nmf_overlap,o5_nmf_svd_cor=o5_nmf_svd_cor))

}

validate_marker_std_output_parameter <- function(data.matrix=data.matrix,tProp=tProp,R1_filter_step1_results_new=R1_filter_step1_results_new,NMF_input1=NMF_input1,immune_cell_uni_table0_GS=immune_cell_uni_table0_GS,cor.cutoff=0.9,comb=comb,
									list4.1=R1_selectedCM_step2_results_Hclust[[1]],list4.2=R1_selectedCM_step2_results_fixedCT[[1]],list4.3=R1_selectedCM_step2_results_HP[[1]] ){

	output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff=0.9)
	output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff=0.9)

	o1_top100_ICTDmarker = output2[["ccc_overlap_order"]]
	o1_top100_cor_gene = output2[["ccc_whole_order"]]
	o2_R1_overlap = output1[["common_mat"]]
	o3_R4_overlap = output2[["common_mat"]]
	o3_R4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]])
	o4_R4_4.1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.1)
	o4_R4_4.2_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.2)
	o4_R4_4.3_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.3)
	o5_nmf_overlap = output2[["common_NMF"]]
	o5_nmf_svd_cor = output2[["cor_top_row"]]

	return(list(o1_top100_ICTDmarker=o1_top100_ICTDmarker,o1_top100_cor_gene=o1_top100_cor_gene,o2_R1_overlap=o2_R1_overlap,
				o3_R4_overlap=o3_R4_overlap,o3_R4_svd_cor=o3_R4_svd_cor,
				o4_R4_4.1_svd_cor=o4_R4_4.1_svd_cor,o4_R4_4.2_svd_cor=o4_R4_4.2_svd_cor,o4_R4_4.3_svd_cor=o4_R4_4.3_svd_cor,
				o5_nmf_overlap=o5_nmf_overlap,o5_nmf_svd_cor=o5_nmf_svd_cor))

}


validate_marker_std_output_parameter_4.4 <- function(data.matrix=data.matrix,tProp=tProp,R1_filter_step1_results_new=R1_filter_step1_results_new,NMF_input1=NMF_input1,immune_cell_uni_table0_GS=immune_cell_uni_table0_GS,cor.cutoff=0.9,comb=comb,
									list4.1=R1_selectedCM_step2_results_Hclust[[1]],list4.2=R1_selectedCM_step2_results_fixedCT[[1]],list4.3=R1_selectedCM_step2_results_HP[[1]] ,
									list4.4=R1_selectedCM_step2_results_HCTES[[1]]){

	output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)
	output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)

	o1_top100_ICTDmarker = output2[["ccc_overlap_order"]]
	o1_top100_cor_gene = output2[["ccc_whole_order"]]
	o2_R1_overlap = output1[["common_mat"]]
	o2_R1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]])
	o3_R4_overlap = output2[["common_mat"]]
	o3_R4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]])
	o4_R4_4.1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.1)
	o4_R4_4.2_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.2)
	o4_R4_4.3_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.3)
	o4_R4_4.4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.4)
	o5_nmf_overlap = output2[["common_NMF"]]
	o5_nmf_svd_cor = output2[["cor_top_row"]]

	return(list(o1_top100_ICTDmarker=o1_top100_ICTDmarker,o1_top100_cor_gene=o1_top100_cor_gene,
				o2_R1_overlap=o2_R1_overlap,o2_R1_svd_cor=o2_R1_svd_cor,
				o3_R4_overlap=o3_R4_overlap,o3_R4_svd_cor=o3_R4_svd_cor,
				o4_R4_4.1_svd_cor=o4_R4_4.1_svd_cor,o4_R4_4.2_svd_cor=o4_R4_4.2_svd_cor,o4_R4_4.3_svd_cor=o4_R4_4.3_svd_cor, o4_R4_4.4_svd_cor=o4_R4_4.4_svd_cor,
				o5_nmf_overlap=o5_nmf_overlap,o5_nmf_svd_cor=o5_nmf_svd_cor))

}

validate_marker_std_output_parameter_4.5 <- function(data.matrix=data.matrix,tProp=tProp,R1_filter_step1_results_new=R1_filter_step1_results_new,NMF_input1=NMF_input1,immune_cell_uni_table0_GS=immune_cell_uni_table0_GS,cor.cutoff=0.9,
									list4.1=R1_selectedCM_step2_results_Hclust[[1]],list4.2=R1_selectedCM_step2_results_fixedCT[[1]],list4.3=R1_selectedCM_step2_results_HP[[1]] ,
									list4.4=R1_selectedCM_step2_results_HCTES[[1]],list4.5=tg_selected_R4_RR){

	output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)
	output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)

	o1_top100_ICTDmarker = output2[["ccc_overlap_order"]]
	o1_top100_cor_gene = output2[["ccc_whole_order"]]
	o2_R1_overlap = output1[["common_mat"]]
	o2_R1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]])
	o3_R4_overlap = output2[["common_mat"]]
	o3_R4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]])
	o4_R4_4.1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.1)
	o4_R4_4.2_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.2)
	o4_R4_4.3_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.3)
	o4_R4_4.4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.4)
	o4_R4_4.5_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.5)
	o5_nmf_overlap = output2[["common_NMF"]]
	o5_nmf_svd_cor = output2[["cor_top_row"]]

	return(list(o1_top100_ICTDmarker=o1_top100_ICTDmarker,o1_top100_cor_gene=o1_top100_cor_gene,
				o2_R1_overlap=o2_R1_overlap,o2_R1_svd_cor=o2_R1_svd_cor,
				o3_R4_overlap=o3_R4_overlap,o3_R4_svd_cor=o3_R4_svd_cor,
				o4_R4_4.1_svd_cor=o4_R4_4.1_svd_cor,o4_R4_4.2_svd_cor=o4_R4_4.2_svd_cor,o4_R4_4.3_svd_cor=o4_R4_4.3_svd_cor, o4_R4_4.4_svd_cor=o4_R4_4.4_svd_cor,o4_R4_4.5_svd_cor=o4_R4_4.5_svd_cor,
				o5_nmf_overlap=o5_nmf_overlap,o5_nmf_svd_cor=o5_nmf_svd_cor))

}

validate_marker_std_output_parameter_4.6_4.7 <- function(data.matrix=data.matrix,tProp=tProp,R1_filter_step1_results_new=R1_filter_step1_results_new,NMF_input1=NMF_input1,immune_cell_uni_table0_GS=immune_cell_uni_table0_GS,cor.cutoff=0.9,
									list4.1=R1_selectedCM_step2_results_Hclust[[1]],list4.2=R1_selectedCM_step2_results_fixedCT[[1]],list4.3=R1_selectedCM_step2_results_HP[[1]] ,
									list4.4=R1_selectedCM_step2_results_HCTES[[1]],list4.5=tg_selected_R4_RR,
									list4.6=CTES3, list4.7=NMF_selected_R1){

	#output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)
	#output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1[[2]],table=immune_cell_uni_table0_GS,cor.cutoff)

output1 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]],NMFmarker=NMF_input1,table=immune_cell_uni_table0_GS,cor.cutoff)	
output2 = validate_marker_coverage_by_proportion(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]],NMFmarker=NMF_input1,table=immune_cell_uni_table0_GS,cor.cutoff)


	o1_top100_ICTDmarker = output2[["ccc_overlap_order"]]
	o1_top100_cor_gene = output2[["ccc_whole_order"]]
	o2_R1_overlap = output1[["common_mat"]]
	o2_R1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[1]])
	o3_R4_overlap = output2[["common_mat"]]
	o3_R4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=R1_filter_step1_results_new[[4]])
	o4_R4_4.1_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.1)
	o4_R4_4.2_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.2)
	o4_R4_4.3_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.3)
	o4_R4_4.4_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.4)
	o4_R4_4.5_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.5)
	o4_R4_4.6_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.6)
	o4_R4_4.7_svd_cor = validata_R4(origiData=data.matrix,truePro=tProp,markList=list4.7)
	o5_nmf_overlap = output2[["common_NMF"]]
	o5_nmf_svd_cor = output2[["cor_top_row"]]

	return(list(o1_top100_ICTDmarker=o1_top100_ICTDmarker,o1_top100_cor_gene=o1_top100_cor_gene,
				o2_R1_overlap=o2_R1_overlap,o2_R1_svd_cor=o2_R1_svd_cor,
				o3_R4_overlap=o3_R4_overlap,o3_R4_svd_cor=o3_R4_svd_cor,
				o4_R4_4.1_svd_cor=o4_R4_4.1_svd_cor,o4_R4_4.2_svd_cor=o4_R4_4.2_svd_cor,o4_R4_4.3_svd_cor=o4_R4_4.3_svd_cor, 
				o4_R4_4.4_svd_cor=o4_R4_4.4_svd_cor,o4_R4_4.5_svd_cor=o4_R4_4.5_svd_cor,
				o4_R4_4.6_svd_cor=o4_R4_4.6_svd_cor,o4_R4_4.7_svd_cor=o4_R4_4.7_svd_cor,
				o5_nmf_overlap=o5_nmf_overlap,o5_nmf_svd_cor=o5_nmf_svd_cor))

}

save_R4_marker <- function(R4_list = R1_filter_step1_results_new[[4]]){
	
	R4_unique = list()
	llll_name = c()
	tt = table(names(R4_list))
	for(i in 1:length(tt)){
		llll_name = c(llll_name, names(tt[i]))
		loca = which(names(R4_list)==names(tt[i]) )
		new_flag = 0
		for(j in 1:length(loca)){
			#print(loca[j])
			if(new_flag==0){
				R4_unique[[i]] = R4_list[[loca[j] ]]
				new_flag = 1
			}else{
				R4_unique[[i]] = c(R4_unique[[i]], R4_list[[loca[j] ]])
			}
		}
	}
	names(R4_unique) = llll_name

	R4_unique = lapply(R4_unique, unique)

	return(R4_unique)
}

#######################################
# test 75688 neuron and B cell 

#data_nnn <- data.matrix[R1_filter_step1_results_new[[4]][[1]],]
#data_bbb <- data.matrix[R1_filter_step1_results_new[[4]][[33]],]
#v_nnn = (svd(data_nnn)$v)[, 1]
#v_bbb = (svd(data_bbb)$v)[, 1]

#cor(v_nnn, v_bbb)
#> R1_filter_step1_results_new[[4]][[1]]
# [1] "LRMP"     "MEF2C"    "LCP1"     "STAP1"    "HLA-DMB"  "TCL1A"    "FCRLA"    "KIAA0922" "HMGB1"    "DLGAP5"  
#> R1_filter_step1_results_new[[4]][[33]]
# [1] "LRMP"     "DLGAP5"   "TCL1A"    "FCRLA"    "KIAA0922" "HMGB1"    "DTX1"     "CDKN3"    "DEPDC1"   "TPX2"   


run_NMF_r2 <- function(NMF_input1=NMF_input1, maxIter=20000){

###########Preprocess the data
P_fix=NMF_input1$S_indi
NMF_indi_all = NMF_input1$NMF_indi_all

data_t=NMF_input1$NMF_data
addP=NMF_input1$NMF_P_pre
NMF_indi_all0<-NMF_indi_all[which(apply(NMF_indi_all,1,sum)<=2),]
NMF_indi_all01<-NMF_indi_all0[,which(apply(NMF_indi_all0,2,sum)>0 | P_fix==0)] #?????0#######if 5 is okay, originally it was 10
NMF_indi_all01<-NMF_indi_all01[which(apply(NMF_indi_all01,1,sum)>0),]
NMF_indi_all=NMF_indi_all01[rownames(NMF_indi_all01)%in%rownames(data_t),]
X1=data_t[match(rownames(NMF_indi_all),rownames(data_t)),]###take only those rows that correspond to rows in NMF_indi_all
K=ncol(NMF_indi_all)
indiS=1-NMF_indi_all
indiS[,which(P_fix==0)]=0*indiS[,which(P_fix==0)]
###########
	
###########Parameter settings
theta=500 ##penalty parameter for constraints on NMF_indi_all
indiS_method="nonprdescent" ##the updating scheme for the structural constraints
iter=maxIter
alpha=beta=gamma=roh=0
#gamma=0.1
#roh=0.1 ##roh has to be equal to 0 in order to implement addP
UM=VM=NULL
qq=1
epslog=8
nPerm=2
initial_U=initial_V=NULL
mscale=0###########1
cnormalize=1
###########

###########Run the constrained qNMF
ttt1=qnmf_indisS_all_revise_addP(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale,cnormalize,addP, P_fix)
###########
X1=ttt1$X1
U=ttt1$U
V=ttt1$V
#sum((X1-(ttt1$U)%*%t(ttt1$V))^2)
#diag(cor(ttt1$U,NMF_indi_all, method="spearman"))
###########Run cross validation of the constrained qNMF
#devs1=cv_qnmf(X1,indiS_method,NMF_indi_all,alpha,beta,gamma,roh,theta,qq,iter,epslog,nPerm,mscale,cnormalize)
#############################################

return(list(X1=X1, U=U, V=V, ttt1=ttt1))

}


plot_cor_bar <- function(m_epic=epic.ccc,m_timer=timer.ccc,m_ictd=ictd.ccc,GSEname=dataSetName,check=overlap_mark){
	library(ggplot2)
	library(reshape2)
	library(gridExtra)

	if(ncol(m_epic)!=ncol(m_timer) | ncol(m_epic)!=ncol(m_ictd) )
		stop("dim of correlation matrix incorrect!chang")

	max_timer = apply(m_timer, 2, max)
	max_epic = apply(m_epic, 2, max)
	#max_ictd = apply(m_ictd, 2, max)
	max_4.1 =apply(check[["o4_R4_4.1_svd_cor"]], 1, max)
	max_4.2 =apply(check[["o4_R4_4.2_svd_cor"]], 1, max)
	max_4.3 =apply(check[["o4_R4_4.3_svd_cor"]], 1, max)
	max_R4 =apply(check[["o3_R4_svd_cor"]], 1, max) 

	aaa <- rbind(max_timer,max_epic,max_4.1,max_4.2,max_4.3,max_R4)
	aaa <- as.data.frame(aaa)
	aaa1 <- cbind(rownames(aaa),aaa)
	colnames(aaa1) <- c("method", colnames(aaa))
	#if(GSEname!="GSE75688simu_diri")  aaa1 <- aaa1[,-5]
	vvlibrary <- unlist(strsplit(GSEname, "_", fixed=T))
	if("GSE72056" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE81861" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE103322" %in% vvlibrary ) aaa1 <- aaa1[,-c(3,11)]

	aaa2 <- melt(aaa1, id='method')
	aaa2$method <- factor(aaa2$method, levels=c("max_timer","max_epic","max_4.1","max_4.2","max_4.3","max_R4") )
	cbPalette <- c("burlywood", "cornsilk1", "skyblue1", "skyblue2", "skyblue3", "skyblue4")

		#barplot(aaa[,1], main=paste(GSEname,colnames(m_epic)[i],sep=""),col = c("lightblue"), names.arg="" )
		mm<- ggplot(aaa2, aes(x = variable, y = value, fill = method))+
			geom_bar(stat = 'identity', position = 'dodge')+
			scale_fill_manual(values=cbPalette)+
			 geom_text(aes(label=round(value,3), y = value-0.13), position = position_dodge(1),angle=90,color="black")+
			 ggtitle(paste(GSEname,"correlation comparison",sep="") )+
			 theme(plot.title = element_text(hjust = 0.5))


	#plot_ff = paste("plot_cor_bar_", GSEname, ".pdf", sep="")
	#pdf(plot_ff,width=10,height=5)
	#print(mm)
	#dev.off()
	return(mm)
}

plot_bar_ES <- function(m_epic=epic.ccc,m_timer=timer.ccc,m_ictd=ictd.ccc,GSEname=dataSetName,check=overlap_mark,pp=pp_list){
	library(ggplot2)
	library(reshape2)
	library(gridExtra)

	if(ncol(m_epic)!=ncol(m_timer) | ncol(m_epic)!=ncol(m_ictd) )
		stop("dim of correlation matrix incorrect!chang")

	max_timer = apply(m_timer, 2, max)
	max_epic = apply(m_epic, 2, max)
	#max_ictd = apply(m_ictd, 2, max)
	max_4.1 =apply(check[["o4_R4_4.1_svd_cor"]], 1, max)
	max_4.2 =apply(check[["o4_R4_4.2_svd_cor"]], 1, max)
	max_4.3 =apply(check[["o4_R4_4.3_svd_cor"]], 1, max)
	max_R4 =apply(check[["o3_R4_svd_cor"]], 1, max) 

	aaa <- rbind(max_timer,max_epic,max_4.1,max_4.2,max_4.3,max_R4)
	aaa <- as.data.frame(aaa)
	aaa1 <- cbind(rownames(aaa),aaa)
	colnames(aaa1) <- c("method", colnames(aaa))
	#if(GSEname!="GSE75688simu_diri")  aaa1 <- aaa1[,-5]
	vvlibrary <- unlist(strsplit(GSEname, "_", fixed=T))
	if("GSE72056" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE81861" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE103322" %in% vvlibrary ) aaa1 <- aaa1[,-c(3,11)]

	aaa2 <- melt(aaa1, id='method')
	aaa2$method <- factor(aaa2$method, levels=c("max_timer","max_epic","max_4.1","max_4.2","max_4.3","max_R4") )
	cbPalette <- c("burlywood", "cornsilk1", "skyblue1", "skyblue2", "skyblue3", "skyblue4")

		#barplot(aaa[,1], main=paste(GSEname,colnames(m_epic)[i],sep=""),col = c("lightblue"), names.arg="" )
		mm<- ggplot(aaa2, aes(x = variable, y = value, fill = method))+
			geom_bar(stat = 'identity', position = 'dodge')+
			scale_fill_manual(values=cbPalette)+
			 geom_text(aes(label=round(value,3), y = value-0.13), position = position_dodge(1),angle=90,color="black")+
			 ggtitle(paste(GSEname,"correlation comparison",sep="") )+
			 theme(plot.title = element_text(hjust = 0.5))


	#plot_ff = paste("plot_cor_bar_", GSEname, ".pdf", sep="")
	#pdf(plot_ff,width=10,height=5)
	#print(mm)
	#dev.off()

	pp[[1]] <- mm

	plot_ff = paste("ES_", dataSetName, ".pdf", sep="")
	pdf(plot_ff,width=15,height=10)
	lay <- rbind(c(1,1,1),
				 c(2,3,4))
	#do.call(grid.arrange, c(pp,nrow=1 ) )
	#do.call(grid.arrange(grobs = pp, layout_matrix = lay) )
	grid.arrange(grobs = pp, layout_matrix = lay)
	dev.off()
	
}

plot_bar_ES_addtest <- function(m_epic=epic.ccc,m_timer=timer.ccc,m_ictd=ictd.ccc,GSEname=dataSetName,check=overlap_mark,pp=pp_list){
	library(ggplot2)
	library(reshape2)
	library(gridExtra)

	if(ncol(m_epic)!=ncol(m_timer) | ncol(m_epic)!=ncol(m_ictd) )
		stop("dim of correlation matrix incorrect!chang")

	max_timer = apply(m_timer, 2, max)
	max_epic = apply(m_epic, 2, max)
	#max_ictd = apply(m_ictd, 2, max)
	max_4.1 =apply(check[["o4_R4_4.1_svd_cor"]], 1, max)
	max_4.2 =apply(check[["o4_R4_4.2_svd_cor"]], 1, max)
	max_4.3 =apply(check[["o4_R4_4.3_svd_cor"]], 1, max)
	max_4.4 =apply(check[["o4_R4_4.4_svd_cor"]], 1, max) 
	max_R4 =apply(check[["o3_R4_svd_cor"]], 1, max) 

	aaa <- rbind(max_timer,max_epic,max_4.1,max_4.2,max_4.3,max_4.4, max_R4)
	aaa <- as.data.frame(aaa)
	aaa1 <- cbind(rownames(aaa),aaa)
	colnames(aaa1) <- c("method", colnames(aaa))
	#if(GSEname!="GSE75688simu_diri")  aaa1 <- aaa1[,-5]
	vvlibrary <- unlist(strsplit(GSEname, "_", fixed=T))
	if("GSE72056" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE81861" %in% vvlibrary) aaa1 <- aaa1[,-5]
	if("GSE103322" %in% vvlibrary ) aaa1 <- aaa1[,-c(3,11)]

	aaa2 <- melt(aaa1, id='method')
	aaa2$method <- factor(aaa2$method, levels=c("max_timer","max_epic","max_4.1","max_4.2","max_4.3","max_4.4","max_R4") )
	cbPalette <- c("burlywood", "cornsilk1", "skyblue1", "skyblue2", "skyblue3", "lightblue","skyblue4")

		#barplot(aaa[,1], main=paste(GSEname,colnames(m_epic)[i],sep=""),col = c("lightblue"), names.arg="" )
		mm<- ggplot(aaa2, aes(x = variable, y = value, fill = method))+
			geom_bar(stat = 'identity', position = 'dodge')+
			scale_fill_manual(values=cbPalette)+
			 geom_text(aes(label=round(value,3), y = value-0.13), position = position_dodge(1),angle=90,color="black")+
			 ggtitle(paste(GSEname,"correlation comparison",sep="") )+
			 theme(plot.title = element_text(hjust = 0.5))


	#plot_ff = paste("plot_cor_bar_", GSEname, ".pdf", sep="")
	#pdf(plot_ff,width=10,height=5)
	#print(mm)
	#dev.off()

	pp[[1]] <- mm

	plot_ff = paste("ES_new", dataSetName, ".pdf", sep="")
	pdf(plot_ff,width=15,height=10)
	lay <- rbind(c(1,1,1),
				 c(2,3,4))
	#do.call(grid.arrange, c(pp,nrow=1 ) )
	#do.call(grid.arrange(grobs = pp, layout_matrix = lay) )
	grid.arrange(grobs = pp, layout_matrix = lay)
	dev.off()
	
}


save_data_try_parameter <- function(){

#try parameter
bulk_data = data.matrix
tProp =tProp
R1_whole = R1_filter_step1_results_new
data_CORS_cancer = data_CORS_cancer
overlap_mark = overlap_mark

#cancer_module = cancer_module
#module_selected = module_selected
MR_IM_result_new_c = MR_IM_result_new_c
list_new_c2 = list_new_c2
R1_filter_step1_results_new=R1_filter_step1_results_new
R1_selectedCM_step2_results_new=R1_selectedCM_step2_results_new
#NMF_input1=NMF_input1

R1_selectedCM_step2_results_Hclust
R1_selectedCM_step2_results_fixedCT
R1_selectedCM_step2_results_HP

f_name=paste(dataSetName,"_parameter.RData",sep="")
save(list=c("bulk_data","tProp","R1_whole","data_CORS_cancer","overlap_mark",
	"cancer_module","module_selected","MR_IM_result_new_c","list_new_c2","R1_filter_step1_results_new","R1_selectedCM_step2_results_new","NMF_input1"), 
	file=f_name)

}



check_gene_expr_in_cell <- function(sc_RNA, Cell_Type, genelist,GSEid){

	print(paste("genelist has total ",length(genelist), "gene",sep=""))
  	print(paste("# of common gene between SC and genelist is", length(intersect(rownames(sc_RNA),genelist)),sep=""))
 # genelist=intersect(rownames(sc_RNA),genelist)
  genelist=intersect(genelist,rownames(sc_RNA))
  #sc_RNA = GSE103322_scRNA
  #Cell_Type = GSE103322_cell_type
  #unique(Cell_Type)
  #table(Cell_Type)

  if(GSEid=="GSE81861") cell_library <- c("Bcell","Fibroblast","Tcell","Endothelial","Epithelial","Macrophage","MastCell","unknown")
  if(GSEid=="GSE75688") cell_library <- c("6","1","0","2","3","4","5")
  if(GSEid=="GSE72056") cell_library <- c("B_cell","CAF","T_cell","Endo","Macro","NK","Malignant","unknown")
  if(GSEid=="GSE103322")cell_library <- c("B_cell","CAF","T_cell","Endothelial","Macrophage","Dendritic","Myofibro","Malignant1","Malignant2")

  data_subset = list()
  for(i in 1:length(cell_library)){
  		
      col_loca = which(Cell_Type==cell_library[[i]] )
      data_subset[[i]] <- sc_RNA[genelist, col_loca, drop=F]
  } 
  names(data_subset) <- cell_library

  data_mean = list()
  for(i in 1:length(cell_library)){
  	
      data_mean[[i]] <- apply(data_subset[[i]],1,mean)
  }
  names(data_mean) <- cell_library

 ttt <- c()
for(i in 1:length(data_mean))
{
	ttt <- cbind(ttt,data_mean[[i]])
}
colnames(ttt) <- names(data_mean)

  return(list(data_subset=data_subset, data_mean_matrix=ttt) )
}

check_gene_expr_in_cell_new <- function(sc_RNA, Cell_Type, genelist,GSEid){

	print(paste("genelist has total ",length(genelist), "gene",sep=""))
  	print(paste("# of common gene between SC and genelist is", length(intersect(rownames(sc_RNA),genelist)),sep=""))
 # genelist=intersect(rownames(sc_RNA),genelist)
  genelist=intersect(genelist,rownames(sc_RNA))
  #sc_RNA = GSE103322_scRNA
  #Cell_Type = GSE103322_cell_type
  #unique(Cell_Type)
  #table(Cell_Type)

  if(GSEid=="GSE81861") cell_library <- c("Bcell","Fibroblast","Tcell","Endothelial","Epithelial","Macrophage","MastCell","unknown")
  if(GSEid=="GSE75688") cell_library <- names(table(Cell_Type))
  if(GSEid=="GSE72056") cell_library <- c("B_cell","CAF","T_cell","Endo","Macro","NK","Malignant","unknown")
  if(GSEid=="GSE103322")cell_library <- c("B_cell","CAF","T_cell","Endothelial","Macrophage","Dendritic","Myofibro","Malignant1","Malignant2")

  data_subset = list()
  for(i in 1:length(cell_library)){
  		
      col_loca = which(Cell_Type==cell_library[[i]] )
      data_subset[[i]] <- sc_RNA[genelist, col_loca, drop=F]
  } 
  names(data_subset) <- cell_library

  data_mean = list()
  for(i in 1:length(cell_library)){
  	
      data_mean[[i]] <- apply(data_subset[[i]],1,mean)
  }
  names(data_mean) <- cell_library

 ttt <- c()
for(i in 1:length(data_mean))
{
	ttt <- cbind(ttt,data_mean[[i]])
}
colnames(ttt) <- names(data_mean)

  return(list(data_subset=data_subset, data_mean_matrix=ttt) )
}


melt_for_boxplot <- function(cell_rmse=list)
{

	yyy = c()
	for(i in 1:length(cell_rmse)){
		len = length(cell_rmse[[i]])
		vvv = matrix(cell_rmse[[i]],len,1)
		ggg = matrix(names(cell_rmse[[i]]),len,1)
		cc = rep(names(cell_rmse)[i],len)
		mmm = matrix(cc,len,1)
		xxx = cbind.data.frame(mmm, ggg, vvv)
		yyy = rbind(yyy, xxx)
	}

	return(yyy)

}

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



group_CIBER_marker <- function(mylist = cell_marker_c)
{

	#length(cell_marker_c) #22
	B_CIBERmarker <- c(mylist[["B.cells.naive"]], mylist[["B.cells.memory"]])
	CD4T_CIBERmark <- c(mylist[["T.cells.CD4.naive"]],mylist[["T.cells.CD4.memory.resting"]],mylist[["T.cells.CD4.memory.activated"]],mylist[["T.cells.follicular.helper"]],mylist[["T.cells.regulatory..Tregs."]],mylist[["T.cells.gamma.delta"]])
	CD8T_CIBERmark <- c(mylist[["T.cells.CD8"]])
	NK_CIBERmark <- c(mylist[["NK.cells.resting"]],mylist[["NK.cells.activated"]])
	Monocyte_CIBERmark <- mylist[["Monocytes"]]
	Macrophage_CIBERmark <- c(mylist[["Macrophages.M0"]],mylist[["Macrophages.M1"]],mylist[["Macrophages.M2"]])
	Dendritic_CIBERmark <- c(mylist[["Dendritic.cells.resting"]],mylist[["Dendritic.cells.activated"]])
	Mast_CIBERmark <- c(mylist[["Mast.cells.resting"]],mylist[["Mast.cells.activated"]])
	Neutrophils_CIBERmark <- mylist[["Neutrophils"]]

	CIBER_mark <- list(B_CIBERmarker,CD4T_CIBERmark,CD8T_CIBERmark,NK_CIBERmark,Monocyte_CIBERmark,Macrophage_CIBERmark,Dendritic_CIBERmark,Mast_CIBERmark,Neutrophils_CIBERmark)
	names(CIBER_mark) <- c("B_CIBERmarker","CD4T_CIBERmark","CD8T_CIBERmark","NK_CIBERmark","Monocyte_CIBERmark","Macrophage_CIBERmark","Dendritic_CIBERmark","Mast_CIBERmark","Neutrophils_CIBERmark")

	return(CIBER_mark)
}


make_venn_matrix <- function(inputlist = marker.4method)
{
	#length(inputlist)
	#names(inputlist)
	all_gene <- c()
	for(i in 1:length(inputlist))
	{
		all_gene <- c(all_gene, inputlist[[i]])
	}
	all_gene <- unique(all_gene)

	mmm <- c()
	aaa <- data.frame(gene = all_gene, aaa=1)
	ictd <- data.frame(gene = inputlist[[1]], ictd=1)
	timer <- data.frame(gene = inputlist[[2]], timer =1)
	epic <- data.frame(gene = inputlist[[3]], epic =1)
	ciber <- data.frame(gene= inputlist[[4]],ciber=1)
	true_mk <- data.frame(gene = inputlist[[5]], true_mk=1)
	ictd_r1 <- data.frame(gene = inputlist[[6]], ictd_r1=1)

#	mmm <- merge(aaa, ictd,by.x="gene", by.y="gene", all=T)
#	mmm <- merge(mmm, timer, by.x = "gene", by.y="gene", all=T)
#	mmm <- merge(mmm, epic, by.x = "gene", by.y="gene", all = T)
#	mmm <- merge(mmm, ciber, by.x = "gene", by.y="gene", all = T)
#	mmm <- merge(mmm, true_mk, by.x = "gene", by.y="gene", all = T)

	mmm <- merge(aaa, ictd, all=T)
	mmm <- merge(mmm, timer, all=T)
	mmm <- merge(mmm, epic, all = T)
	mmm <- merge(mmm, ciber,  all = T)
	mmm <- merge(mmm, true_mk, all = T)
	mmm <- merge(mmm, ictd_r1, all = T)

	mmm[is.na(mmm)] <- 0

	return(mmm)
}

extend_marker_from_R1 <- function(cell_marker_i, R1_filter_step1_results_new)
{

	new_cell_mark <- list()
	for(i in 1:length(cell_marker_i))
	{
		cell_mk_lib <- c()
		for(j in 1:length(cell_marker_i[[i]]))
		{
			vv <- cell_marker_i[[i]][[j]]

			cell_mk_lib <- c(cell_mk_lib , R1_filter_step1_results_new[[1]][[vv]])
			cell_mk_lib <- unique(cell_mk_lib)

		}
		new_cell_mark[[i]] <- cell_mk_lib


	}
	names(new_cell_mark) <- names(cell_marker_i)

	return(new_cell_mark)
}


