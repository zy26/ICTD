##library(clue)
#'@importFrom clue solve_LSAP
pMatrix.min <- function(A, B) {##########not useful anymore
####This is to shuffle the rows of A so that A could be as close to B as possible;
    # finds the permutation P of A such that ||PA - B|| is minimum in Frobenius  norm
    # Uses the linear-sum assignment problem (LSAP) solver in the "clue" package

    # Returns P%*%A and the permutation vector `pvec' such that
    # A[pvec, ] is the permutation of A closest to B
    n <- nrow(A)
    D <- matrix(NA, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            D[j, i] <- (sum((B[j, ] - A[i, ])^2))
        } }
    vec <- c(clue::solve_LSAP(D))
    list(A=A[vec,], pvec=vec)
}

#'@importFrom clue solve_LSAP
re_order<-function(U,NMF_indi_all){
	if(ncol(U)>1 & ncol(NMF_indi_all)>0){
		C=cor(NMF_indi_all,U, method="spearman")
		C=C-min(C,na.rm=TRUE)+0.1
		C[is.na(C)]=0
		if(ncol(U)<ncol(NMF_indi_all)){
			C=C[1:ncol(U),1:ncol(U)]
		}else{
			C=cbind(C,matrix(0, nrow=nrow(C),ncol=(ncol(U)-ncol(NMF_indi_all))))
		}
		vec <- c(clue::solve_LSAP(C,maximum = TRUE))
	}else{
		vec=1:ncol(U)
	}
	list(U=U[,vec,drop=FALSE], pvec=vec)
}
nmf1<-function(X){
        S=matrix(0,nrow=ncol(X), ncol=ncol(X))
        for(i in 1:nrow(X)){
                x=X[i,,drop=FALSE]
                S=S+t(x)%*%x
        }
        dec  = abs(eigen(S)$vectors[,1])
        names(dec)=colnames(X)
        return(dec)
}


#'@importFrom NMF basis
calc.initial.UV_shuff<-function(X1,K,NMF_indi_all){
#####function to initialize the U and V matrix for qNMF, and while an indicator matrix, NMF_indi_all, exists, would also shuffle the columns of U and V based on its similarity with the indicator matrix.
#        aaa=kmeans(X1,centers=K,nstart=5)
#        initial_U=t(sapply(aaa$cluster, function(x)(c(1:K)==x)*1))
#        aaa=kmeans(X1,centers=K,nstart=5)
#        initial_V=t(aaa$centers)

#        aaa=kmeans(t(X1),centers=K,nstart=5)
#        initial_U=t(aaa$centers)
#        aaa=kmeans(X1,centers=K,nstart=5)
#        initial_V=t(aaa$centers)

#library(NMF)
	if(K>1){
		res = NMF::nmf(X1, K)
		initial_U=NMF::basis(res)
		initial_V=t(NMF::coef(res))
	}else{
		res=nmf1(X1)
		initial_V=matrix(res,ncol=1)
		initial_U=compute.V(t(X1),(initial_V))
	}

        A=t(initial_U)
        B=t(NMF_indi_all) ###B should be the matrix where 1 means genes to keep, while 0 means genes to not keep;
        #X <- pMatrix.min(A,B)
        X <- re_order(initial_U,NMF_indi_all)
        initial_U=initial_U[,X$pvec,drop=FALSE]
        initial_V=initial_V[,X$pvec,drop=FALSE]
	for(j in 1:ncol(initial_U)){
		initial_V[,j]=initial_V[,j]*sqrt(sum(initial_U[,j]^2))
		initial_U[,j]=initial_U[,j]/sqrt(sum(initial_U[,j]^2))
	}
        aaa=list(initial_U, initial_V)
        names(aaa)=c("initial_U", "initial_V")
        return(aaa)
}
#'@importFrom NMF basis
calc.initial.FSG_shuff<-function(X1,K1,K2,NMF_indi_all){
#####function to initialize the F and S matrix for tNMF, and while an indicator matrix, NMF_indi_all, exists, would also shuffle the columns of F and G based on its similarity with the indicator matrix.
#	aaa=kmeans(X1,centers=K1,nstart=5)
#	initial_F=t(sapply(aaa$cluster, function(x)(c(1:K1)==x)*1))
#	aaa=kmeans(t(X1),centers=K1,nstart=5)
#	initial_F=t(aaa$centers)
#	aaa=kmeans(X1,centers=K2,nstart=5)
#	initial_G=t(aaa$centers)
#library(NMF)
	if(K>1){
		res = NMF::nmf(X1, K1)
		initial_F=NMF::basis(res)
		initial_G=t(NMF::coef(res))
	}else{
		res=nmf1(X1)
		initial_G=matrix(res,ncol=1)
		initial_F=compute.V(t(X1),(initial_G))
	}

        A=t(initial_F)
        B=t(NMF_indi_all) ###B should be the matrix where 1 means genes to keep, while 0 means genes to not keep;
        #X <- pMatrix.min(A,B)
        X <- re_order(initial_F,NMF_indi_all)
        initial_F=initial_F[,X$pvec,drop=FALSE]
        initial_G=initial_G[,X$pvec,drop=FALSE]

	initial_S=(t(initial_F)%*%X1%*%initial_G)
	initial_S=diag(ncol(initial_F))
	aaa=list(initial_F, initial_S,initial_G)
	names(aaa)=c("initial_F", "initial_S", "initial_G")
	return(aaa)
}






initialize_mats<-function(input,K,K_common,ddim,NMF_indi_all){
	if(class(input)=="list"){
		K_vec=K
		X1_list=input
		X1_all=NULL
		for(i in 1:length(X1_list)){
		        X1=X1_list[[i]]
		        X1_all=cbind(X1_all, X1)
		}
		if(ddim=="qnmf"){
			initial_U_list=initial_V_list=vector("list", length(X1_list))
			NMF_indi_all_list=vector("list", length(X1_list))
		        aaaa=calc.initial.UV_shuff(X1_all,K_common,NMF_indi_all)
		        U_common=aaaa[[1]]
			print("The common U")
			print(head(U_common))
		        V_tmp=aaaa[[2]]
		        iiii = cbind(c(1,cumsum(sapply(X1_list, ncol))[-length(X1_list)]+1), cumsum(sapply(X1_list, ncol)))
		        for(i in 1:length(X1_list)){
				NMF_indi_all_list[[i]]=NMF_indi_all
		                X1=X1_list[[i]]
				if(K_common<K_vec[i]){
					X1_remain=X1-(U_common%*%t(V_tmp[iiii[i,1]:iiii[i,2],]))
					X1_remain=(X1_remain>0)*X1_remain
					bbbb=calc.initial.UV_shuff(X1_remain,K_vec[i]-K_common,NMF_indi_all[,-c(1:K_common),drop=FALSE])
					U_spec=bbbb[[1]]
				}else{
					U_spec=NULL
				}
				initial_U=cbind(U_common,U_spec)
				for(j in 1:ncol(initial_U)){
					initial_U[,j]=initial_U[,j]/sqrt(sum((initial_U[,j])^2))
				}
				initial_V=compute.V(X1, initial_U)
				#print("The initial U,V")
				#print(head(initial_U))
				#print(head(initial_V))
		                initial_U_list[[i]]=initial_U
		                initial_V_list[[i]]=initial_V
		        }
		        print("step2 ini")

			ini_list=list(initial_U_list,initial_V_list)
			names(ini_list)=c("initial_U_list","initial_V_list")
		}else{
			initial_F_list=initial_S_list=initial_G_list=vector("list", length(X1_list))
			NMF_indi_all_list=vector("list", length(X1_list))
		        aaaa=calc.initial.UV_shuff(X1_all,K_common,NMF_indi_all)
		        F_common=aaaa[[1]]
		        G_tmp=aaaa[[2]]
		        iiii = cbind(c(1,cumsum(sapply(X1_list, ncol))[-length(X1_list)]+1), cumsum(sapply(X1_list, ncol)))
			for(i in 1:length(X1_list)){
				NMF_indi_all_list[[i]]=NMF_indi_all
			        X1=X1_list[[i]]
				if(K_common<K_vec[i]){
					X1_remain=X1-(F_common%*%t(G_tmp[iiii[i,1]:iiii[i,2],]))
					X1_remain=(X1_remain>0)*X1_remain
					bbbb=calc.initial.UV_shuff(X1_remain,K_vec[i]-K_common,NMF_indi_all[,-c(1:K_common),drop=FALSE])
					F_spec=bbbb[[1]]
				}else{
					F_spec=NULL
				}
				initial_F=cbind(F_common,F_spec)
				initial_G=compute.V(X1, initial_F)
				initial_S=diag(ncol(initial_F))
			        initial_F_list[[i]]=initial_F
			        initial_S_list[[i]]=initial_S
			        initial_G_list[[i]]=initial_G
			
			}
			ini_list=list(initial_F_list,initial_S_list,initial_G_list)
			names(ini_list)=c("initial_F_list","initial_S_list","initial_G_list")
		}

	}else{
		X1=input
		if(ddim=="qnmf"){
			UV=calc.initial.UV_shuff(X1,K, NMF_indi_all)
			initial_U=UV[[1]]
			initial_V=compute.V(X1, initial_U)
			#print("test V")
			#print(head(initial_V))
			for(j in 1:ncol(initial_U)){
				initial_V[,j]=initial_V[,j]*sqrt(sum((initial_U[,j])^2))
				initial_U[,j]=initial_U[,j]/sqrt(sum((initial_U[,j])^2))
			}
			#print("The initial U,V")
			#print(head(initial_U))
			#print(head(initial_V))	#bug: V has wrong colname
			
			ini_list=list(initial_U,initial_V)
			names(ini_list)=c("initial_U","initial_V")
		}else{
			FSG=calc.initial.FSG_shuff(X1,K,K, NMF_indi_all)
			initial_F=FSG[[1]]
			initial_S=FSG[[2]]
			initial_G=FSG[[3]]
			#######

			ini_list=list(initial_F,initial_S,initial_G)
			names(ini_list)=c("initial_F","initial_S","initial_G")
		}
	}
	return(ini_list)
}
initialize_mats1<-function(input,K,ddim,NMF_indi_all){
	if(class(input)=="list"){
		X1_list=input
		if(ddim=="qnmf"){
			initial_U_list=initial_V_list=vector("list", length(X1_list))
			for(i in 1:length(X1_list)){
				NMF_indi_all_list[[i]]=NMF_indi_all
			        X1=X1_list[[i]]
			
			        UV=calc.initial.UV_shuff(X1,K, NMF_indi_all)
				initial_U=UV[[1]]; initial_V=UV[[2]];
			
			        initial_U_list[[i]]=initial_U
			        initial_V_list[[i]]=initial_V
			
			}
			ini_list=list(initial_U_list,initial_V_list)
			names(ini_list)=c("initial_U_list","initial_V_list")
		}else{
			initial_F_list=initial_S_list=initial_G_list=vector("list", length(X1_list))
			NMF_indi_all_list=vector("list", length(X1_list))
			for(i in 1:length(X1_list)){
				NMF_indi_all_list[[i]]=NMF_indi_all
			        X1=X1_list[[i]]
			        FSG=calc.initial.FSG_shuff(X1,K,K, NMF_indi_all)
				initial_F=FSG[[1]]; initial_S=FSG[[2]]; initial_G=FSG[[3]]
			#	initial_F=initial_F*NMF_indi_all_list[[i]]
			#	initial_F[initial_F<0.01]=0.01
			
			        initial_F_list[[i]]=initial_F
			        initial_S_list[[i]]=initial_S
			        initial_G_list[[i]]=initial_G
			
			}
			ini_list=list(initial_F_list,initial_S_list,initial_G_list)
			names(ini_list)=c("initial_F_list","initial_S_list","initial_G_list")
		}

	}else{
		X1=input
		if(ddim=="qnmf"){
			UV=calc.initial.UV_shuff(X1,K, NMF_indi_all)
			initial_U=UV[[1]]
			initial_V=UV[[2]]
			
			ini_list=list(initial_U,initial_V)
			names(ini_list)=c("initial_U","initial_V")
		}else{
			FSG=calc.initial.FSG_shuff(X1,K,K, NMF_indi_all)
			initial_F=FSG[[1]]
			initial_S=FSG[[2]]
			initial_G=FSG[[3]]
			#######

			ini_list=list(initial_F,initial_S,initial_G)
			names(ini_list)=c("initial_F","initial_S","initial_G")
		}
	}
	return(ini_list)
}


initialize_GL<-function(input,K,ddim){ ##initialize the Graph constraints
	if(class(input)=="list"){
		X1_list=input
		if(ddim=="qnmf"){
			UM_list=DU_list=VM_list=DV_list=NULL
			for(i in 1:length(X1_list)){
			        X1=X1_list[[i]]
			
				G=nrow(X1); N=ncol(X1)
				UM=matrix(0,nrow=G, ncol=G)
				DU=diag(rowSums(UM))
				VM=matrix(0,nrow=N, ncol=N)
				DV=diag(rowSums(VM))
			        DV_list[[i]]=DV
			        UM_list[[i]]=UM
			        DU_list[[i]]=DU
			        VM_list[[i]]=VM
			}
			ini_list=list(UM_list,DU_list,VM_list, DV_list)
			names(GL_list)=c("UM_list","DU_list","VM_list","DV_list")
		}else{
			FM_list=DF_list=GM_list=DG_list=NULL
			for(i in 1:length(X1_list)){
			        X1=X1_list[[i]]
				
				G=nrow(X1); N=ncol(X1)
				FM=matrix(0,nrow=G, ncol=G)
				DF=diag(rowSums(FM))
				GM=matrix(0,nrow=N, ncol=N)
				DG=diag(rowSums(GM))
			        DG_list[[i]]=DG
			        FM_list[[i]]=FM
			        DF_list[[i]]=DF
			        GM_list[[i]]=GM
			}
			GL_list=list(FM_list,DF_list,GM_list,DG_list)
			names(GL_list)=c("FM_list","DF_list","GM_list","DG_list")
		}

	}else{
		X1=input
		G=nrow(X1); N=ncol(X1)
		if(ddim=="qnmf"){
			UM=matrix(0,nrow=G, ncol=G)
			DU=diag(rowSums(UM))
			VM=matrix(0,nrow=N, ncol=N)
			DV=diag(rowSums(VM))
			GL_list=list(UM, DU, VM, DV)
			names(GL_list)=c("UM","DU","VM","DV")
		}else{
			FM=matrix(0,nrow=G, ncol=G)
			DF=diag(rowSums(FM))
			GM=matrix(0,nrow=N, ncol=N)
			DG=diag(rowSums(GM))
			GL_list=list(FM,DF,GM,DG)
			names(GL_list)=c("FM","DF","GM","DG")
		}
	}
	return(GL_list)
}
