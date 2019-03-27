#2018-May-05


#source("./estimateScore_0.1.r")
#source("./coexpression_WGCNA_0.2.r")




#circumstance:
# a. cancer, normal, purity,
# b. cancer, purity,
# c. cancer, normal
# d. cancer


find_cancer_module <- function(str,Fdata_t, Fdata_n="NO", purity="NO")
{
	#backup data
	Fdata_t0 <- Fdata_t
	Fdata_n0 <- Fdata_n

	#save different flag, track information
	track_information <- c("Without normal data, no DEG step")
	length1 <- 0.1
	length2 <- 0.1

	# prepare 1:
	# judge the data type, if max value > 30, then take log(x+1)
	if(max(Fdata_t0) > 30){
		Fdata_t0 <- log(Fdata_t0 + 1)
	}


	# prepare 2:
	# replace ENS ID with gene symbol
	#tg_ids<-TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!="") ,7]#only protein coding
	#Fdata_t0 <- Fdata_t[tg_ids,]
	#rownames(Fdata_t0) <- TCGA_ensem_annotation[which(TCGA_ensem_annotation[,3]=="protein_coding" & TCGA_ensem_annotation[,4]!=""),4]
	#rownames(Fdata_t0) <- make.names(rownames(Fdata_t0),unique=T)
	write.table(Fdata_t0, "./data.txt",sep="\t",quote=F)

	#==========================================================================
	#  (i) based on cancer mixture data (normal data, purity data are optional),
	#     We get genes which highly correlation with purity.
	#    If no purity data, using ESTIMATE get ourself purity.
	#    If provide normal data, we can narrow the range of gene by picking up up-regulated genes.
	#	Do not know how to use CCLE and cancer_common data???
	#==========================================================================

	# 1. judge purity flag, if do not have purity information, calculate it.
	if(purity == "NO"){
		# using ESTIMATE to give purity
		filterCommonGenes(input.f="./data.txt",output.f="data.gct",id="GeneSymbol")
		myscore <- estimateScore_0.1("data.gct","data_score.gct",platform="affymetrix")
		colnames(myscore) <- gsub(".","-",colnames(myscore),fixed=T)
		my_purity <- myscore[c(-1,-2,-3),]	#just left the cancer purity score
	}else{
		# the user already give the purity data,
		# we clean the purity data and transfor it into myscore for doing correlation.
		stop("please provide purity information")
	}

	# 2. then, correlated with above purity in order to give gene name
	# parameter: threshold of correlation
	#myscoreV <- as.vector(myscore)
	my_correlation <- cor(t(Fdata_t0),t(my_purity))
	cor_sort <- sort(my_correlation[,1],decreasing = TRUE)#[1:1000]
	correlation_threshold <- 0
	corre_gene <- cor_sort[which(cor_sort > correlation_threshold)]
	gene_number_after_correlated <- length(corre_gene)
	data_m <- Fdata_t0[names(corre_gene),]

	# 3. use CCLE  (CCLE_expressed_genes)
	#length(intersect(rownames(data_t0),CCLE_gene_symbol )  )


	# 4. judge normal flag
	if( Fdata_n != "NO" ){
		# using normal data to narrow the range which just upregulated with cancer
		Fdata_n0 <- Fdata_n[tg_ids,]

		DEG_stat <- wilcox_test_all(Fdata_t0,Fdata_n0)
		up_regulated_genes <- rownames(DEG_stat)[which((DEG_stat[,4]<0.05)&(DEG_stat[,2]>0))] #FDR less than 0.05 and sign is positive
		down_regulated_gene <-  rownames(DEG_stat)[which((DEG_stat[,4]<0.05)&(DEG_stat[,2]<0))]
		whole_set <- rownames(data_m)
		up_and_down <- union(up_regulated_genes, down_regulated_gene)
		not_change_gene <- setdiff(whole_set, up_and_down)
		not_down <- union(up_regulated_genes, not_change_gene)

		length1 <- length(intersect(up_regulated_genes,CCLE_gene_symbol))
		length2 <- length(intersect(not_down,CCLE_gene_symbol))

		if(length(intersect(up_regulated_genes,CCLE_gene_symbol))>2000){
			gene_name_adj <- intersect(up_regulated_genes,CCLE_gene_symbol)	#use up regulated gene
			track_information <- c("use up-regulated to narrow gene range")

		}else if(length(intersect(not_down,CCLE_gene_symbol)) >2000){
			gene_name_adj <- intersect(not_down,CCLE_gene_symbol) #use not down regulated gene
			track_information <- c("use not-down-regulated to narrow gene range")
		}else{
			gene_name_adj <-  rownames(data_m)		#do not use Deg because common gene is not enough
			track_information <- c("Do not use DEG to narrow (not enough gene),even provide normal data")
		}



		gene_tmp <- intersect(rownames(data_m), gene_name_adj)
		data_m <- data_m[gene_tmp,] # gene number decrease from 60000 to 7448
	}


	#====================================================================
	#
	# (ii) Using WGCNA to get the co-expressive module
	#
	#====================================================================

	cancerType <- "str"
	cancer_module <- c()
	row_number <- nrow(data_m)
	col_number <- ncol(data_m)
	#track_information <- c(track_information, input_to_WGCNA_gene_number=row_number, input_to_WGCNA_gene_number=col_number)
	cancer_module <- coexpression_WGCNA_0.2(data_m, cancerType)

	#==================
	# return list
	#==================

	other_information = list(gene_number_after_correlated=gene_number_after_correlated,
			track_information=track_information,
			NumberOf_up_intersect_CCLE=length1, NumberOf_not_down_intersect_CCLE=length2,
			input_to_WGCNA_gene_number=row_number, input_to_WGCNA_sample_number=col_number)

	obj <- list(cancer_module=cancer_module, ESTIMATE_purity=myscore, other_information=other_information)

	return(obj)	#cancer module (this is a list)

}




