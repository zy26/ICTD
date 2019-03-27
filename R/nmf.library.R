tnmf_indisS_all_revise<-function(X1,initial_F,initial_S, initial_G,NMF_indi_all,indiS_method,FM,GM,alpha,beta,gamma,roh,theta,qq, iter, epslog, mscale, cnormalize){
##########This is for tri matrix factorization, indiS is an indicator matrix equal to 1-NMF_indi_all; indiS_method has two options, either project descent, or gradient descent; theta is the penalty on similarity of F and NMF_indi_all matrix; qq is the value controling the updating scheme; 
##X1
##initial_F
##initial_S
##initial_G
##NMF_indi_all
##indiS_method
##FM
##GM
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
##mscale
#library(pracma)
	K_common=NULL
	if(mscale==1){
		X1=t(apply(X1, MARGIN=1, function(x)x/median(x)))
	}
	indiS=1-NMF_indi_all
	if(is.matrix(initial_F) & is.matrix(initial_S) & is.matrix(initial_G)){##both initial F and G are given
	}else{
		inis.list.tnmf=initialize_mats(X1,ncol(NMF_indi_all),K_common,"tnmf",NMF_indi_all)
		initial_F=inis.list.tnmf[["initial_F"]]
		initial_S=inis.list.tnmf[["initial_S"]]
		initial_G=inis.list.tnmf[["initial_G"]]
	}

	GL.list.tnmf = initialize_GL(X1,K,"tnmf")
	if(is.matrix(FM)){
		DF=diag(rowSums(FM))
	}else{
		FM=GL.list.tnmf[["FM"]]
		DF=diag(rowSums(FM))
	}
	if(is.matrix(GM)){
		DG=diag(rowSums(GM))
	}else{
		GM=GL.list.tnmf[["GM"]]
		DG=GL.list.tnmf[["DG"]]
	}

	F=initial_F; G=initial_G; S=initial_S;
####Check the initial residuals
	Xpred=F%*%S%*%t(G)
	print("This is the first iteration")
	print(paste("The fitting residual is", sum((Xpred-X1)^2)))
####Check the initial residuals
	objs=NULL; objs.parts=NULL
	flag=0
	rrr=0
	eps=10^(-epslog) 
#	update_F.list=vector("list",6)
	#iter=6
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_F = F;old_S=S; old_G = G;
		if(gamma==0){
			updateF="HALS-col"
		}else{
			updateF="multiplicative"
		}
		if(updateF=="multiplicative"){
			#optimize F by fixing G;
			a_F = (X1%*%G%*%t(S) + alpha*initial_F + gamma*FM%*%F);
			b_F = (F%*%S%*%t(G)%*%G%*%t(S)+alpha*F+ gamma*DF%*%F);
	                if(indiS_method=="prdescent"){
				b_F[b_F<eps]=eps
	                        update_F = (a_F/b_F)^qq;
	                        update_F = update_F*NMF_indi_all
	                }else{
	                        b_F=b_F+theta*indiS
				b_F[b_F<eps]=eps
	                        update_F = (a_F/b_F)^qq;
	                }
			F[F<eps]=eps #This value is very important! Because in multiplicative update, if you get a zero, it is not going to be updated;
			#update_F[update_F<eps]=eps
			F  = F * update_F;
			#Normalization
		}
		if(updateF=="HALS-col"){
			V=G%*%t(S)
			for(i in 1:ncol(F)){
				tmp.mat=NULL
				for(j in 1:ncol(G)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*F[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_F=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_F[,i,drop=FALSE]
				b_F=(sum(V[,i]*V[,i])+alpha)
				b_F[b_F<eps]=eps
		                if(indiS_method=="prdescent"){
					update_F=a_F/b_F
		                        update_F = update_F*(NMF_indi_all[,i]) ###force those entries that are zero in indiS matrix to be zero
		                }else{
                			a_F=a_F-theta*indiS[,i]
					update_F=a_F/b_F
		                }
				update_F = update_F*(update_F>0)
				if(cnormalize==1){
					update_F[update_F<eps]=eps
					S[i,]=S[i,]*sqrt(sum(update_F^2))##############!!!!!!!!!!!???????
					update_F=update_F/sqrt(sum(update_F^2))
				}
				#update_F[update_F<eps]=eps
				F[,i]=update_F
			}
		}
		if(roh==0){
			updateG="HALS-col"
		}else{
			updateG="multiplicative"
		}
		if(updateG=="multiplicative"){
			#optimize G by fixing F;
			a_G = (t(X1)%*%F%*%S +beta*initial_G + roh*GM%*%G);
			b_G = (G%*%t(S)%*%t(F)%*%F%*%S +beta*G+ roh*DG%*%G);
			b_G[b_G<eps]=eps
			update_G = (a_G/b_G)^qq;
			#update_G = (a_G/b_G);
			#update_G[update_G<eps]=eps
			G[G<eps]=eps
			G  = G * update_G;
			#Normalization
		}
		if(updateG=="HALS-col"){
			U=F%*%S
			for(i in 1:ncol(G)){
				tmp.mat=NULL
				for(j in 1:ncol(F)){
					if(j!=i){
						tmp=sum(U[,j]*U[,i])*G[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_G=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_G[,i,drop=FALSE]
				b_G=(sum(U[,i]*U[,i])+beta)
				b_G[b_G<eps]=eps
				update_G=a_G/b_G
				update_G = update_G*(update_G>0)
				#update_G[update_G<eps]=eps
				G[,i]=update_G
			}
		}
		#optimize S by fixing F and G;
		a_S = t(F)%*%X1%*%G;
		b_S = t(F)%*%F%*%S%*%t(G)%*%G;
		#b_S = t(F)%*%F%*%S%*%t(G)%*%G; ##correction based on Chris Ding paper
		b_S[b_S<eps]=eps
		#update_S = sqrt(a_S/b_S);
		update_S = (a_S/b_S);
		S[S<eps]=eps
		S = S * update_S;

		#obj_tmp=sum((X1-F%*%S%*%G)^2)+alpha*sum((F-initial_F)^2)+beta*sum((G-initial_G)^2)+
		s1=sum((X1-F%*%S%*%t(G))^2)
		s21=sum((F-initial_F)^2); s22=sum((G-initial_G)^2)
		s2=s21*alpha+s22*beta
		s31=sum(diag(t(F)%*%(DF-FM)%*%F)); s32=sum(diag(t(G)%*%(DG-GM)%*%G))
		s3=s31*gamma+s32*roh
		s4=2*theta*sum(diag(t(F)%*%indiS))
		obj_tmp=s1+s2+s3+s4
		objs=c(objs,obj_tmp)
		objs.parts=rbind(objs.parts, matrix(c(s1,s21,s22,s31,s32,s4),nrow=1))

		if(is.finite(obj_tmp)){
			if(rrr>5){
				if(abs(objs[length(objs)]-objs[length(objs)-1])<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		
		Xpred=F%*%S%*%t(G)
		print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	}
	rmse=sum((Xpred-X1)^2)
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", obj_tmp))
	#print(paste("The parameters are",paste(c(alpha, beta, gamma, roh),collapse=" ")))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter,  are",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
	ttt=list(initial_F,initial_S,initial_G,F,S,G,objs,objs.parts,rmse)
	names(ttt)=c("initial_F","initial_S","initial_G","F","S","G","objs","objs.parts","rmse")
	return(ttt)

}


tnmf_indisS_sub_revise<-function(X1,initial_F,initial_S, initial_G,NMF_indi_all,indiS_method,FM,GM,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common, cnormalize){
##########This is simialr to tnmf_indisS_all_revise, only that updateF only allows for "HALS-cols", and only the K-K_common columns will be updated in F; and it always only iterate for once
##X1
##initial_F
##initial_S
##initial_G
##NMF_indi_al
##indiS_method
##FM
##DF
##GM
##DG
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
##K_common
#library(pracma)
	F=initial_F;
	G=initial_G;
	S=initial_S;
	eps=10^(-epslog) 
#	update_F.list=vector("list",6)
	#iter=6
	DF=diag(rowSums(FM))
	DG=diag(rowSums(GM))
	old_F = F;old_S=S; old_G = G;
	if(gamma==0){
		updateF="HALS-col"
	}else{
		updateF="multiplicative"
	}
	indiS=1-NMF_indi_all
	if(updateF=="HALS-col"){ #In MTL, we could only update in a column-wise fashion
		V=G%*%t(S)
		for(i in 1:ncol(F)){
			if(i>K_common){
				tmp.mat=NULL
				for(j in 1:ncol(G)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*F[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_F=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_F[,i,drop=FALSE]
				b_F=(sum(V[,i]*V[,i])+alpha)
				b_F[b_F<eps]=eps
		                if(indiS_method=="prdescent"){
					update_F=a_F/b_F
		                        update_F = update_F*(NMF_indi_all[,i]) ###force those entries that are zero in indiS matrix to be zero
		                }else{
                			a_F=a_F-theta*indiS[,i]
					update_F=a_F/b_F
		                }
				update_F = update_F*(update_F>0)
				
				if(cnormalize==1){
					update_F[update_F<eps]=eps
					S[i,]=S[i,]*sqrt(sum(update_F^2))##############!!!!!!!!!!!???????
					update_F=update_F/sqrt(sum(update_F^2))
				}
				#update_F[update_F<eps]=eps
				F[,i]=update_F
			}
		}
	}
	if(roh==0){
		updateG="HALS-col"
	}else{
		updateG="multiplicative"
	}
	if(updateG=="multiplicative"){
		#optimize G by fixing F;
		a_G = (t(X1)%*%F%*%S +beta*initial_G + roh*GM%*%G);
		b_G = (G%*%t(S)%*%t(F)%*%F%*%S +beta*G+ roh*DG%*%G);
		b_G[b_G<eps]=eps
		update_G = (a_G/b_G)^qq;
		#update_G = (a_G/b_G);
	#	update_G[update_G<eps]=eps
		G[G<eps]=eps
		G  = G * update_G;
		#Normalization
	}
	if(updateG=="HALS-col"){
		U=F%*%S
		for(i in 1:ncol(G)){
			tmp.mat=NULL
			for(j in 1:ncol(F)){
				if(j!=i){
					tmp=sum(U[,j]*U[,i])*G[,j]
					tmp.mat=cbind(tmp.mat, tmp)
				}
			}
			a_G=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_G[,i,drop=FALSE]
			b_G=(sum(U[,i]*U[,i])+beta)
			b_G[b_G<eps]=eps
			update_G=a_G/b_G
			update_G = update_G*(update_G>0)
#			update_G[update_G<eps]=eps
			G[,i]=update_G
		}
	}
	#optimize S by fixing F and G;
	a_S = t(F)%*%X1%*%G;
	b_S = t(F)%*%F%*%S%*%t(G)%*%G;
	#b_S = t(F)%*%F%*%S%*%t(G)%*%G; ##correction based on Chris Ding paper
	b_S[b_S<eps]=eps
	update_S = sqrt(a_S/b_S);
	#update_S = (a_S/b_S);
	S[S<eps]=eps
	S = S * update_S;
	#obj_tmp=sum((X1-F%*%S%*%G)^2)+alpha*sum((F-initial_F)^2)+beta*sum((G-initial_G)^2)+

	ttt=list(F,S,G)
	names(ttt)=c("F","S","G")
	return(ttt)
}



tnmf_indisS_all_MTL<-function(X1_list,initial_F_list,initial_S_list, initial_G_list,NMF_indi_all_list,indiS_method,FM_list,GM_list,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common,mscale,cnormalize){
#####This is the main function for MTL, and it is similar to tnmf_indisS_all_revise, except that K_common columns of F are updated differently using all the tasks
##X1_list
##initial_F_list
##initial_S_list
##initial_G_list
##NMF_indi_all_list
##indiS_method
##FM_list
##DF_list
##GM_list
##DG_list
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##updateF
##updateG
##epslog
##K_common
	if(is.list(initial_F_list)){
	}else{
		initial_F_list=vector("list", length(X1_list))
	}
	if(is.list(initial_S_list)){
	}else{
		initial_S_list=vector("list", length(X1_list))
	}
	if(is.list(initial_G_list)){
	}else{
		initial_G_list=vector("list", length(X1_list))
	}
	if(is.list(FM_list)){
	}else{
		FM_list=vector("list", length(X1_list))
	}
	if(is.list(GM_list)){
	}else{
		GM_list=vector("list", length(X1_list))
	}

####add to tMTL
	X1_all=NULL
	for(i in 1:length(X1_list)){
		X1=X1_list[[i]]
		X1_all=cbind(X1_all, X1)
	}
	X1_med=apply(X1_all, MARGIN=1, median,na.rm=TRUE)
####add to tMTL
        for(i in 1:length(X1_list)){
                NMF_indi_all = NMF_indi_all_list[[i]]
                X1=X1_list[[i]]
	        if(mscale==1){
        	        X1=(apply(X1, MARGIN=2, function(x)x/X1_med))
			X1_list[[i]]=X1
	        }
		initial_F=initial_F_list[[i]]
		initial_S=initial_S_list[[i]]
		initial_G=initial_G_list[[i]]
		if(is.matrix(initial_F) &is.matrix(initial_S) & is.matrix(initial_G)){##both initial F and G are given
		}else{
			inis.list.tnmf=initialize_mats(X1,ncol(NMF_indi_all),"tnmf",NMF_indi_all)
			initial_F=inis.list.tnmf[["initial_F"]]
			initial_S=inis.list.tnmf[["initial_S"]]
			initial_G=inis.list.tnmf[["initial_G"]]
		}
####add to tMTL
		if(cnormalize==1){
			for(j in 1:K_common){
				initial_G[,j]=initial_G[,j]*sum((initial_F[,j])^2)
				initial_F[,j]=initial_F[,j]/sum((initial_F[,j])^2)
			}
		}
####add to tMTL
                initial_F_list[[i]]=initial_F
                initial_S_list[[i]]=initial_S
                initial_G_list[[i]]=initial_G
		GL.list.tnmf = initialize_GL(X1,K,"tnmf")
		FM=FM_list[[i]]
		if(is.matrix(FM)){
		}else{
			FM=GL.list.tnmf[["FM"]]
		}
		GM=GM_list[[i]]
		if(is.matrix(GM)){
		}else{
			GM=GL.list.tnmf[["GM"]]
		}
		FM_list[[i]]=FM
		GM_list[[i]]=GM
        }

	if(theta==0){
		for(i in 2:length(X1_list)){
			initial_F=initial_F_list[[i]]
			initial_G=initial_G_list[[i]]
		        X <- re_order(initial_F,initial_F_list[[1]])
		        initial_F=initial_F[,X$pvec]
        		initial_G=initial_G[,X$pvec]
                	initial_F_list[[i]]=initial_F
	                initial_G_list[[i]]=initial_G
		}
	}
####add to tMTL
#	for(i in 1:length(X1_list)){
#		U_common=U_common+initial_U_list[[i]][,1:K_common]*ncol(X1_list[[i]])
#	}
#	U_common=U_common/sum(sapply(X1_list, ncol))
#	for(i in 1:length(X1_list)){
#		X1=X1_list[[i]]
#		initial_U=initial_U_list[[i]]
#		initial_V=initial_V_list[[i]]
#		initial_U[,1:K_common]=U_common
#		initial_V=compute.V(X1, initial_U)
#        	initial_U_list[[i]]=initial_U
#                initial_V_list[[i]]=initial_V
#	}
	X1_all=NULL
	for(i in 1:length(X1_list)){
		X1=X1_list[[i]]
		X1_all=cbind(X1_all, X1)
	}
	aaaa=calc.initial.FSG_shuff(X1_all,K_common,NMF_indi_all)
	F_common=aaaa[[1]]
	G_tmp=aaaa[[2]]
	print("step1 ini")
	iiii = cbind(c(1,cumsum(sapply(X1_list, ncol))[-length(X1_list)]+1), cumsum(sapply(X1_list, ncol)))
	for(i in 1:length(X1_list)){
		X1=X1_list[[i]]
		initial_F=initial_F_list[[i]]
		initial_G=initial_G_list[[i]]
		initial_F[,1:K_common]=F_common
		initial_G=G_tmp[iiii[i,1]:iiii[i,2],]
        	initial_F_list[[i]]=initial_F
                initial_G_list[[i]]=initial_G
	}
	print("step2 ini")

####add to tMTL

	flag=0
	rrr=0
	eps=10^(-epslog) 
#	update_F.list=vector("list",6)
	#iter=6
	F_list=initial_F_list;S_list=initial_S_list; G_list=initial_G_list;
	objs_mat=NULL
####Check the initial residuals
	fit_res=NULL
	for(k in 1:length(X1_list)){
		X1=X1_list[[k]]
		F=F_list[[k]]
		S=S_list[[k]]
		G=G_list[[k]]
		Xpred=F%*%S%*%t(G)
		fit_res=c(fit_res,sum((Xpred-X1)^2))
	}
	fit_res=c(fit_res, sum(fit_res))
	print("This is the first iteration")
	print(paste("The fitting residual is", paste(fit_res,collapse=", ")))
####Check the initial residuals
	
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_F_list=F_list; old_S_list=S_list; old_G_list=G_list

		for(i in 1:K_common){
			a_F_combine=matrix(0,nrow=nrow(F),ncol=1)
			b_F_combine=0
			for(k in 1:length(X1_list)){
				initial_F=initial_F_list[[k]]
				initial_S=initial_S_list[[k]]
				initial_G=initial_G_list[[k]]
				X1=X1_list[[k]]
				F=F_list[[k]]
				S=S_list[[k]]
				G=G_list[[k]]
				NMF_indi_all=NMF_indi_all_list[[k]]
				indiS=1-NMF_indi_all
				FM=FM_list[[k]]
				DF=diag(rowSums(FM))
				GM=GM_list[[k]]
				DG=diag(rowSums(GM))
				tmp.mat=NULL
				V=G%*%t(S)
				for(j1 in 1:ncol(G)){
					if(j1!=i){
						tmp=sum(V[,j1]*V[,i])*F[,j1]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_F=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_F[,i,drop=FALSE]
				b_F=(sum(V[,i]*V[,i])+alpha)
				a_F_combine=a_F_combine+a_F
				b_F_combine=b_F_combine+b_F
			}
			b_F_combine[b_F_combine<eps]=eps
	                if(indiS_method=="prdescent"){
				update_F=a_F_combine/b_F_combine
	                        update_F = update_F*(NMF_indi_all[,i]) ###force those entries that are zero in indiS matrix to be zero
	                }else{
        			a_F_combine=a_F_combine-theta*indiS[,i]
				update_F=a_F_combine/b_F_combine
	                }
			update_F = update_F*(update_F>0)
			#update_F[update_F<eps]=eps
			if(cnormalize==1){
				update_F[update_F<eps]=eps
				S_list[[k]][i,]=S_list[[k]][i,]*sqrt(sum(update_F^2))##############!!!!!!!!!!!???????
				update_F=update_F/sqrt(sum(update_F^2))
			}
			for(k in 1:length(X1_list)){
				tmp.mat=NULL
				F_list[[k]][,i]=update_F
			}
		}
		for(j in 1:length(X1_list)){
			initial_F=initial_F_list[[j]]
			initial_S=initial_S_list[[j]]
			initial_G=initial_G_list[[j]]
			X1=X1_list[[j]]
			F=F_list[[j]]
			S=S_list[[j]]
			G=G_list[[j]]
			NMF_indi_all=NMF_indi_all_list[[j]]
			indiS=1-NMF_indi_all
			FM=FM_list[[j]]
			DF=diag(rowSums(FM))
			GM=GM_list[[j]]
			DG=diag(rowSums(GM))
			ttt=tnmf_indisS_sub_revise(X1,F,S, G,NMF_indi_all,indiS_method,FM,GM,alpha,beta,gamma,roh,theta,qq, 1, epslog, K_common,cnormalize)
#tnmf_indisS_sub_revise<-function(X1,initial_F,initial_S, initial_G,NMF_indi_all,indiS_method,FM,DF,GM,DG,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common){
			F=ttt[[1]];S=ttt[[2]];G=ttt[[3]]
			F_list[[j]]=F;S_list[[j]]=S;G_list[[j]]=G

		}

		#obj_tmp=sum((X1-F%*%S%*%G)^2)+alpha*sum((F-initial_F)^2)+beta*sum((G-initial_G)^2)+
		objs=NULL
		for(k in 1:length(X1_list)){
			X1=X1_list[[k]]
			F=F_list[[k]]
			S=S_list[[k]]
			G=G_list[[k]]
			initial_F=initial_F_list[[k]]
			initial_S=initial_S_list[[k]]
			initial_G=initial_G_list[[k]]
			NMF_indi_all=NMF_indi_all_list[[k]]
			indiS=1-NMF_indi_all
			FM=FM_list[[k]]
			DF=diag(rowSums(FM))
			GM=GM_list[[k]]
			DG=diag(rowSums(GM))
			s1=sum((X1-F%*%S%*%t(G))^2)
			s21=sum((F-initial_F)^2); s22=sum((G-initial_G)^2)
			s2=s21*alpha+s22*beta
			s31=sum(diag(t(F)%*%(DF-FM)%*%F)); s32=sum(diag(t(G)%*%(DG-GM)%*%G))
			s3=s31*gamma+s32*roh
			s4=2*theta*sum(diag(t(F)%*%indiS))
			obj_tmp=s1+s2+s3+s4
			objs=c(objs,obj_tmp)
	
		}
		objs_mat=rbind(objs_mat, c(objs,sum(objs)))
		fit_res=NULL
		for(k in 1:length(X1_list)){
			X1=X1_list[[k]]
			F=F_list[[k]]
			S=S_list[[k]]
			G=G_list[[k]]
			Xpred=F%*%S%*%t(G)
			fit_res=c(fit_res,sum((Xpred-X1)^2))
		}
		fit_res=c(fit_res, sum(fit_res))
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(max(abs(objs_mat[nrow(objs_mat),]-objs_mat[nrow(objs_mat)-1,]))<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		print(paste("The fitting residual is", paste(fit_res,collapse=", ")))
	}
	rmse=fit_res[length(fit_res)]
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", objs_mat[nrow(objs_mat),ncol(objs_mat)]))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
	#print(paste("The parameters are",paste(c(alpha, beta, gamma, roh),collapse=" ")))
	#print(paste("The residual is ",paste(c(res1,res2,res3),collapse=" ")))
	#print(paste("The parts objectives are ",paste(c(s1, s21,s22,s31,s32),collapse=" ")))
	ttt=list(initial_F_list,initial_S_list,initial_G_list,F_list,S_list,G_list,objs_mat,rmse)
	names(ttt)=c("initial_F_list","initial_S_list","initial_G_list","F_list","S_list","G_list","objs_mat","rmse")
	return(ttt)
}

compute.V<-function(X,U){
        r=ncol(U)
        res.mat=X*0
        coff.mat=res.mat[1:r,,drop=FALSE]
        for(k in 1:ncol(X)){
                y=X[,k,drop=FALSE]
                fit=nnls(U, y)
                alpha=coef(fit)
                coff.mat[,k]=alpha
        }
	V=t(coff.mat)
        return(V)
}

qnmf_indisS_all_revise<-function(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq, iter, epslog,mscale, cnormalize){
##########This is for qua matrix factorization, indiS is an indicator matrix equal to 1-NMF_indi_all; indiS_method has two options, either project descent, or gradient descent; theta is the penalty on similarity of U and NMF_indi_all matrix; qq is the value controling the updating scheme; 
##X1: original matrix 
##initial_U
##initial_V
##NMF_indi_all
##indiS_method
##UM
##DU
##VM
##DV
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
#library(pracma)
	K_common=NULL
        if(mscale==1){
                X1=t(apply(X1, MARGIN=1, function(x)x/median(x)))
        }
	indiS=1-NMF_indi_all
	if(is.matrix(initial_U) & is.matrix(initial_V)){##both initial U and V are given
	}else{
		if(is.matrix(initial_U)){ ##only initial U is given
			initial_V=compute.V(X1, initial_U)
		}else if(is.matrix(initial_V)){ ##only initial V is given
			initial_U=t(compute.V(t(X1),initial_V))
		}else{ ##none of initial U or V is given
			inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all),K_common,"qnmf",NMF_indi_all)
			initial_U=inis.list.qnmf[["initial_U"]]
			initial_V=inis.list.qnmf[["initial_V"]]
		}
	}
	GL.list.qnmf = initialize_GL(X1,K,"qnmf")
	if(is.matrix(UM)){
		DU=diag(rowSums(UM))
	}else{
		UM=GL.list.qnmf[["UM"]]
		DU=diag(rowSums(UM))
	}
	if(is.matrix(VM)){
		DV=diag(rowSums(VM))
	}else{
		VM=GL.list.qnmf[["VM"]]
		DV=diag(rowSums(VM))
	}

	U=initial_U; V=initial_V;
####Check the initial residuals
	Xpred=U%*%t(V)
	print("This is the first iteration")
	print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	print("This is the second iteration")
####Check the initial residuals
	objs=NULL; objs.parts=NULL
	flag=0
	rrr=0
	eps=10^(-epslog) 
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_U = U;old_V = V;
		if(gamma>0){
			updateU="multiplicative"
		}else{
			updateU="HALS-col"
		}
		if(updateU=="multiplicative"){
			#optimize U by fixing V;
			a_U = (X1%*%V + alpha*initial_U + gamma*UM%*%U);
			b_U = (U%*%t(V)%*%V+alpha*U+ gamma*DU%*%U);
	                if(indiS_method=="prdescent"){
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                        update_U = update_U*NMF_indi_all
	                }else{
	                        b_U=b_U+theta*indiS
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                }
			U[U<eps]=eps #This value is very important! Because in multiplicative update, if you get a zero, it won't be updated at any further step;
			#update_U[update_U<eps]=eps
			U  = U * update_U;
			#Normalization
		}
		if(updateU=="HALS-col"){
		#	U[U<eps]=eps
			for(i in 1:ncol(U)){
				tmp.mat=NULL
				for(j in 1:ncol(V)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*U[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
		                if(indiS_method=="prdescent"){
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in NMF_indi_all matrix to be zero
		                }else{
                			a_U=a_U-theta*indiS[,i]
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                }
				update_U = update_U*(update_U>0)
				if(cnormalize==1){ ####
					update_U[update_U<eps]=eps ####
					V[,i]=V[,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
					update_U=update_U/sqrt(sum(update_U^2)) ####
				} ###
				#update_U[update_U<eps]=eps
				U[,i]=update_U
			}
		}
		if(roh>0){
			updateV="multiplicative"
		}else{
			updateV="HALS-col"
		}
		if(updateV=="multiplicative"){
			#optimize V by fixing U;
			a_V = (t(X1)%*%U +beta*initial_V + roh*VM%*%V);
			b_V = (V%*%t(U)%*%U +beta*V+ roh*DV%*%V);
			b_V[b_V<eps]=eps
			update_V = (a_V/b_V)^qq;
			#update_V[update_V<eps]=eps
			#update_V = (a_V/b_V);
			V[V<eps]=eps
			V  = V * update_V;
			#Normalization
		}
		if(updateV=="HALS-col"){
		#	V[V<eps]=eps
			for(i in 1:ncol(V)){
				tmp.mat=NULL
				for(j in 1:ncol(U)){
					if(j!=i){
						tmp=sum(U[,j]*U[,i])*V[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_V=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_V[,i,drop=FALSE]
				b_V=(sum(U[,i]*U[,i])+beta)
				b_V[b_V<eps]=eps
				update_V=a_V/b_V
				update_V = update_V*(update_V>0)
#			print(update_V)
				update_V[update_V<eps]=eps
				V[,i]=update_V
			}
		}
		
		s1=sum((X1-U%*%t(V))^2)
		s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
		s2=s21*alpha+s22*beta
		s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
		s3=s31*gamma+s32*roh
		s4=2*theta*sum(diag(t(U)%*%indiS))
		obj_tmp=s1+s2+s3+s4

		objs=c(objs,obj_tmp)
		#objs.parts=rbind(objs.parts,s1, s21,s22,s31,s32,s4)
		objs.parts=rbind(objs.parts, matrix(c(s1,s21,s22,s31,s32,s4),nrow=1))
		#print(res)
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(abs(objs[length(objs)]-objs[length(objs)-1])<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		
		Xpred=U%*%t(V)
		print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	}
	rmse=sum((Xpred-X1)^2)
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", obj_tmp))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are:",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
#	print(paste("The residual is ",res,collapse=" "))
	#print(paste("The residual is ",paste(c(res1,res2),collapse=" ")))
	#print(paste("The parts objectives are ",paste(zzz[-1],collapse=" ")))
	ttt=list(initial_U,initial_V,U,V,objs,objs.parts,rmse,X1)
	names(ttt)=c("initial_U","initial_V","U","V","objs","objs.parts","rmse","X1")
	return(ttt)
}


qnmf_indisS_sub_revise<-function(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common,cnormalize){
##########This is simialr to qnmf_indisS_all_revise, only that updateU only allows for "HALS-cols", and only the K-K_common columns will be updated in U; and it always only iterate for once
##X1
##initial_U
##initial_V
##NMF_indi_all
##indiS_method
##UM
##DU
##VM
##DV
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
##K_common
#library(pracma)
	U=initial_U; V=initial_V;
	eps=10^(-epslog)
	DU=diag(rowSums(UM)); DV=diag(rowSums(VM)) 
	old_U = U;old_V = V;
	if(gamma>0){
		updateU="multiplicative"
	}else{
		updateU="HALS-col"
	}
	indiS=1-NMF_indi_all
	if(updateU=="HALS-col"){
		for(i in 1:ncol(U)){
			if(i>K_common){
				tmp.mat=NULL
				for(j in 1:ncol(V)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*U[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
		                if(indiS_method=="prdescent"){
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in indiS matrix to be zero
		                }else{
					if(ncol(indiS)>=i){
	                			a_U=a_U-theta*indiS[,i]
					}
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                }
				update_U = update_U*(update_U>0)
				if(cnormalize==1){
					update_U[update_U<eps]=eps
					V[,i]=V[,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
					update_U=update_U/sqrt(sum(update_U^2))
				}
			#	update_U[update_U<eps]=eps
				U[,i]=update_U
			}
		}
	}
	if(roh>0){
		updateV="multiplicative"
	}else{
		updateV="HALS-col"
	}
	if(updateV=="multiplicative"){
		#optimize V by fixing U;
		a_V = (t(X1)%*%U +beta*initial_V + roh*VM%*%V);
		b_V = (V%*%t(U)%*%U +beta*V+ roh*DV%*%V);
		b_V[b_V<eps]=eps
		update_V = (a_V/b_V)^qq;
		#update_V = (a_V/b_V);
		V[V<eps]=eps
		V  = V * update_V;
		#Normalization
	}
	if(updateV=="HALS-col"){
		for(i in 1:ncol(V)){
			tmp.mat=NULL
			for(j in 1:ncol(U)){
				if(j!=i){
					tmp=sum(U[,j]*U[,i])*V[,j]
					tmp.mat=cbind(tmp.mat, tmp)
				}
			}
			a_V=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_V[,i,drop=FALSE]
			b_V=(sum(U[,i]*U[,i])+beta)
			b_V[b_V<eps]=eps
			update_V=a_V/b_V
			update_V = update_V*(update_V>0)
			#update_V[update_V<eps]=eps
			V[,i]=update_V
		}
	}
	ttt=list(U,V)
	names(ttt)=c("U","V")
	return(ttt)
}



#################need to fix this folloiwng
qnmf_indisS_all_MTL<-function(X1_list,K_vec, initial_U_list, initial_V_list,NMF_indi_all_list,indiS_method,UM_list,VM_list,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common,mscale,cnormalize){###Here, MTL could not work when gamma is larger than 0
#####This is the main function for MTL, and it is similar to qnmf_indisS_all_revise, except that K_common columns of U are updated differently using all the tasks
##X1_list
##initial_U_list
##initial_V_list
##NMF_indi_all_list
##indiS_method
##UM_list
##DU_list
##VM_list
##DV_list
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
##K_common
	flag=0
	rrr=0
	eps=10^(-epslog) 
#	update_U.list=vector("list",6)
	#iter=6

	if(is.list(initial_U_list)){
	}else{
		initial_U_list=vector("list", length(X1_list))
	}
	if(is.list(initial_V_list)){
	}else{
		initial_V_list=vector("list", length(X1_list))
	}
	if(is.list(UM_list)){
	}else{
		UM_list=vector("list", length(X1_list))
	}
	if(is.list(VM_list)){
	}else{
		VM_list=vector("list", length(X1_list))
	}
####add to tMTL
	X1_all=NULL
	for(i in 1:length(X1_list)){
		X1=X1_list[[i]]
		X1_all=cbind(X1_all, X1)
	}
	X1_med=apply(X1_all, MARGIN=1, median,na.rm=TRUE)
####add to tMTL
	    
	if(mscale==1){
        	for(i in 1:length(X1_list)){
                	X1=X1_list[[i]]
        	        X1=(apply(X1, MARGIN=2, function(x)x/X1_med))
			X1_list[[i]]=X1
	        }
	}
	ini_list=initialize_mats(X1_list,K_vec,K_common,"qnmf",NMF_indi_all)
	initial_U_list=ini_list[["initial_U_list"]]
	initial_V_list=ini_list[["initial_V_list"]]
        for(i in 1:length(X1_list)){
		X1=X1_list[[i]]
		initial_V=initial_V_list[[i]]
		initial_U=initial_U_list[[i]]
                NMF_indi_all = NMF_indi_all_list[[i]]
		if(cnormalize==1){
			for(j in 1:K_vec[i]){
				initial_V[,j]=initial_V[,j]*sum((initial_U[,j])^2)
				initial_U[,j]=initial_U[,j]/sum((initial_U[,j])^2)
			}
		}
####add to tMTL
		GL.list.qnmf = initialize_GL(X1,K_vec[i],"qnmf")
		UM=UM_list[[i]]
		if(is.matrix(UM)){
		}else{
			UM=GL.list.qnmf[["UM"]]
		}
		VM=VM_list[[i]]
		if(is.matrix(VM)){
		}else{
			VM=GL.list.qnmf[["VM"]]
		}
		UM_list[[i]]=UM
		VM_list[[i]]=VM
        }


	U_list=initial_U_list; V_list=initial_V_list;
####Check the initial residuals
	fit_res=NULL
	for(k in 1:length(X1_list)){
		X1=X1_list[[k]]
		U=U_list[[k]]
		V=V_list[[k]]
		Xpred=U%*%t(V)
		fit_res=c(fit_res,sum((Xpred-X1)^2))
	}
	fit_res=c(fit_res, sum(fit_res))
	print("This is the first iteration")
	print(paste("The fitting residual is", paste(fit_res,collapse=", ")))
####Check the initial residuals
	
	objs_mat=NULL
	
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_U_list=U_list; old_V_list=V_list
		for(i in 1:K_common){
			a_U_combine=0
			b_U_combine=0
			for(k in 1:length(X1_list)){

				initial_U=initial_U_list[[k]]
				initial_V=initial_V_list[[k]]
				X1=X1_list[[k]]
				U=U_list[[k]]
				V=V_list[[k]]
				NMF_indi_all=NMF_indi_all_list[[k]]
				indiS=1-NMF_indi_all
				UM=UM_list[[k]]
				DU=diag(rowSums(UM))
				VM=VM_list[[k]]
				DV=diag(rowSums(VM))

				tmp.mat=NULL
				for(j1 in 1:ncol(V)){
					if(j1!=i){
						tmp=sum(V[,j1]*V[,i])*U[,j1]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
				b_U=(sum(V[,i]*V[,i])+alpha)
				a_U_combine=a_U_combine+a_U
				b_U_combine=b_U_combine+b_U
			}
			b_U_combine[b_U_combine<eps]=eps
	                if(indiS_method=="prdescent"){
				update_U=a_U_combine/b_U_combine
	                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in NMF_indi_all matrix to be zero
	                }else{
        			a_U_combine=a_U_combine-theta*indiS[,i]
				update_U=a_U_combine/b_U_combine
	                }
			update_U = update_U*(update_U>0)
			if(cnormalize==1){
				for(k in 1:length(X1_list)){
					update_U[update_U<eps]=eps
					V_list[[k]][,i]=V_list[[k]][,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
				}
				update_U=update_U/sqrt(sum(update_U^2))
			}
			#update_U[update_U<eps]=eps
			for(k in 1:length(X1_list)){
				U_list[[k]][,i]=update_U
			}
		}
		for(j in 1:length(X1_list)){
			initial_U=initial_U_list[[j]]
			initial_V=initial_V_list[[j]]
			X1=X1_list[[j]]
			U=U_list[[j]]
			V=V_list[[j]]
			NMF_indi_all=NMF_indi_all_list[[j]]
			indiS=1-NMF_indi_all
			UM=UM_list[[j]]
			DU=diag(rowSums(UM))
			VM=VM_list[[j]]
			DV=diag(rowSums(VM))
			ttt=qnmf_indisS_sub_revise(X1,U, V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq, 1, epslog,K_common,cnormalize)##########
#qnmf_indisS_sub_revise<-function(X1,initial_U,initial_V,indiS,indiS_method,UM,DU,VM,DV,alpha,beta,gamma,roh,theta,qq, iter, epslog, K_common){###########
			U=ttt[[1]];V=ttt[[2]]############
			U_list[[j]]=U;V_list[[j]]=V##########

		}

		#obj_tmp=sum((X1-U%*%V)^2)+alpha*sum((U-initial_U)^2)+beta*sum((V-initial_V)^2)+
		objs=NULL
		for(k in 1:length(X1_list)){
			X1=X1_list[[k]]
			U=U_list[[k]]
			V=V_list[[k]]
			initial_U=initial_U_list[[k]]
			initial_V=initial_V_list[[k]]
			NMF_indi_all=NMF_indi_all_list[[k]]
			indiS=1-NMF_indi_all
			UM=UM_list[[k]]
			DU=diag(rowSums(UM))
			VM=VM_list[[k]]
			DV=diag(rowSums(VM))
			s1=sum((X1-U%*%t(V))^2); s21=sum((U-initial_U)^2)
			s22=sum((V-initial_V)^2)
			s2=s21*alpha+s22*beta
			s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V))
			s3=s31*gamma+s32*roh
			s4=2*theta*sum(diag(t(U)%*%indiS))
			obj_tmp=s1+s2+s3+s4
			objs=c(objs,obj_tmp)
	
		}
		objs_mat=rbind(objs_mat, c(objs,sum(objs)))
		fit_res=NULL
		for(k in 1:length(X1_list)){
			X1=X1_list[[k]]
			U=U_list[[k]]
			V=V_list[[k]]
			Xpred=U%*%t(V)
			fit_res=c(fit_res,sum((Xpred-X1)^2))
		}
		fit_res=c(fit_res, sum(fit_res))
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(max(abs(objs_mat[nrow(objs_mat),]-objs_mat[nrow(objs_mat)-1,]))<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		print(paste("The fitting residual is", paste(fit_res,collapse=", ")))
	}
	rmse=fit_res[length(fit_res)]
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", objs_mat[nrow(objs_mat),ncol(objs_mat)]))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
	#print(paste("The parameters are",paste(c(alpha, beta, gamma, roh),collapse=" ")))
	#print(paste("The residual is ",paste(c(res1,res2,res3),collapse=" ")))
	#print(paste("The parts objectives are ",paste(c(s1, s21,s22,s31,s32),collapse=" ")))
	ttt=list(initial_U_list,initial_V_list,U_list,V_list,objs_mat,rmse)
	names(ttt)=c("initial_U_list","initial_V_list","U_list","V_list","objs_mat","rmse")
	return(ttt)
}

##library(nnls)
#'@importFrom nnls nnls
cv_qnmf<-function(X,indiS_method,NMF_indi_all,alpha,beta,gamma,roh, theta,qq,iter,epslog, nPerm,mscale){
##X
##NMF_indi_all
##theta
##eta
##qq
##iter
##nPerm
##########
        if(mscale==1){
                X=t(apply(X, MARGIN=1, function(x)x/median(x)))
        }
	K=ncol(NMF_indi_all)
	indiS=1-NMF_indi_all
#########
	devs=NULL
        for(rrr in 1:nPerm){
                test_ID = sample(1:ncol(X),0.8*ncol(X))
                X_training=X[,test_ID]
		UM=VM=NULL
		initial_U=NULL
		initial_V=NULL
                X_test=X[,-test_ID]
                ttt=qnmf_indisS_all_revise(X_training,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale,cnormalize)###Run the normal NMF
#qnmf_indisS_all_revise<-function(X,initial_U,initial_V,NMF_indi_all,indiS_method,UM,DU,VM,DV,alpha,beta,gamma,roh,theta,qq, iter, epslog){
		if(is.infinite(ttt[["objs"]][length(ttt[["objs"]])]) | is.nan(ttt[["objs"]][length(ttt[["objs"]])])){
			save(X_training, X_test, file="X_nan.RData")
		}else{
	                U_pred = ttt$U
	                V_pred=NULL
	                for(j in 1:ncol(X_test)){
	                        fit = nnls::nnls(U_pred, X_test[,j])
	                        coffs=coef(fit)
	                        V_pred=rbind(V_pred, coffs)
	                }
	                devs = rbind(devs, c(sd(X_test),quantile(X_test),sum((U_pred%*%t(V_pred)-X_test)^2)))
		}
        }
	colnames(devs)=c("sd","0%","25%", "50%", "75%", "100%","diffs")
        return(devs)
}



#################add true proportion to examin the tendency of theta
#################add fixed rows of P
qnmf_indisS_all_revise_addP_old<-function(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq, iter, epslog,mscale, cnormalize, addP, P_fix){
##########This is for qua matrix factorization, indiS is an indicator matrix equal to 1-NMF_indi_all; indiS_method has two options, either project descent, or gradient descent; theta is the penalty on similarity of U and NMF_indi_all matrix; qq is the value controling the updating scheme; 
##X1: original matrix 
##initial_U
##initial_V
##NMF_indi_all
##indiS_method
##UM
##DU
##VM
##DV
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
#library(pracma)
	K_common=NULL
	K=ncol(NMF_indi_all)
        if(mscale==1){
                X1=t(apply(X1, MARGIN=1, function(x)x/median(x)))
        }
	indiS=1-NMF_indi_all
	indiS[,which(P_fix==0)]=0*indiS[,which(P_fix==0)]
	if(is.matrix(initial_U) & is.matrix(initial_V)){##both initial U and V are given
	}else{
		if(is.matrix(initial_U)){ ##only initial U is given
			initial_V=compute.V(X1, initial_U)
		}else if(is.matrix(initial_V)){ ##only initial V is given
			initial_U=(compute.V(t(X1),initial_V))
		}else{ ##none of initial U or V is given
			MMM=1###MMM=2 is the worst, could pick between 1 and 2
			if(MMM==1){ ##gives better fitting
				NMF_indi_all1=NMF_indi_all
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),ncol(NMF_indi_all),"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=inis.list.qnmf[["initial_V"]]
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}else if(MMM==2){ ##gives better similarity to NMF_indi_all
				NMF_indi_all1=NMF_indi_all[,which(P_fix==1)]
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),K_common,"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=cbind(compute.V(X1, initial_U),t(addP[which(P_fix==0),]))
				initial_U=(compute.V(t(X1),initial_V))
			}else{
				initial_U=(compute.V(t(X1),t(addP[which(P_fix==0),])))
				initial_U=cbind(NMF_indi_all[,which(P_fix==1)],initial_U)
				initial_V=compute.V(X1, initial_U)
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}
		}
	}
	colnames(initial_U)=colnames(NMF_indi_all)
	colnames(initial_V)=colnames(NMF_indi_all)	#chang
	print(dim(initial_V))
	GL.list.qnmf = initialize_GL(X1,K,"qnmf")
	if(is.matrix(UM)){
		DU=diag(rowSums(UM))
	}else{
		UM=GL.list.qnmf[["UM"]]
		DU=diag(rowSums(UM))
	}
	if(is.matrix(VM)){
		DV=diag(rowSums(VM))
	}else{
		VM=GL.list.qnmf[["VM"]]
		DV=diag(rowSums(VM))
	}

	U=initial_U; V=initial_V;
####Check the initial residuals
	Xpred=U%*%t(V)
	print("This is the first iteration")
	print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	print("This is the second iteration")
####Check the initial residuals
	objs=NULL; objs.parts=NULL
	flag=0
	rrr=0
	eps=10^(-epslog) 
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_U = U;old_V = V;
		if(gamma>0){
			updateU="multiplicative"
		}else{
			updateU="HALS-col"
		}
		if(updateU=="multiplicative"){
			#optimize U by fixing V;
			a_U = (X1%*%V + alpha*initial_U + gamma*UM%*%U);
			b_U = (U%*%t(V)%*%V+alpha*U+ gamma*DU%*%U);
	                if(indiS_method=="prdescent"){
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                        update_U = update_U*NMF_indi_all
	                }else{
	                        b_U=b_U+theta*indiS
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                }
			U[U<eps]=eps #This value is very important! Because in multiplicative update, if you get a zero, it won't be updated at any further step;
			#update_U[update_U<eps]=eps
			U  = U * update_U;
			#Normalization
		}
		if(updateU=="HALS-col"){
		#	U[U<eps]=eps
			for(i in 1:ncol(U)){
				tmp.mat=NULL
				for(j in 1:ncol(V)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*U[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
		                if(indiS_method=="prdescent"){
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in NMF_indi_all matrix to be zero
		                }else{
                			a_U=a_U-theta*indiS[,i]
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                }
				update_U = update_U*(update_U>0)
				if(cnormalize==1){ ####
				#if(cnormalize==1 & P_fix[i]==1){ ####
					update_U[update_U<eps]=eps ####
					V[,i]=V[,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
					update_U=update_U/sqrt(sum(update_U^2)) ####
				} ###
				#update_U[update_U<eps]=eps
				U[,i]=update_U
			}
		}
		if(roh>0){
			updateV="multiplicative"
		}else{
			updateV="HALS-col"
		}
		if(updateV=="multiplicative"){
			#optimize V by fixing U;
			a_V = (t(X1)%*%U +beta*initial_V + roh*VM%*%V);
			b_V = (V%*%t(U)%*%U +beta*V+ roh*DV%*%V);
			b_V[b_V<eps]=eps
			update_V = (a_V/b_V)^qq;
			#update_V[update_V<eps]=eps
			#update_V = (a_V/b_V);
			V[V<eps]=eps
			V  = V * update_V;
			#Normalization
		}
		if(updateV=="HALS-col"){
		#	V[V<eps]=eps
			for(i in 1:ncol(V)){
				#if(P_fix[i]==1){
					tmp.mat=NULL
					for(j in 1:ncol(U)){
						if(j!=i){
							tmp=sum(U[,j]*U[,i])*V[,j]
							tmp.mat=cbind(tmp.mat, tmp)
						}
					}
					a_V=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_V[,i,drop=FALSE]
					b_V=(sum(U[,i]*U[,i])+beta)
					b_V[b_V<eps]=eps
					update_V=a_V/b_V
					update_V = update_V*(update_V>0)
	#			print(update_V)
					update_V[update_V<eps]=eps
					V[,i]=update_V
				#}
			}
		}
		
		s1=sum((X1-U%*%t(V))^2)
		s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
		s2=s21*alpha+s22*beta
		s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
		s3=s31*gamma+s32*roh
		s4=2*theta*sum(diag(t(U)%*%indiS))
		obj_tmp=s1+s2+s3+s4

		objs=c(objs,obj_tmp)
		#objs.parts=rbind(objs.parts,s1, s21,s22,s31,s32,s4)
		objs.parts=rbind(objs.parts, matrix(c(s1,s21,s22,s31,s32,s4),nrow=1))
		#print(res)
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(abs(objs[length(objs)]-objs[length(objs)-1])<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		
		Xpred=U%*%t(V)
		print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	}
	rmse=sum((Xpred-X1)^2)
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", obj_tmp))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are:",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
#	print(paste("The residual is ",res,collapse=" "))
	#print(paste("The residual is ",paste(c(res1,res2),collapse=" ")))
	#print(paste("The parts objectives are ",paste(zzz[-1],collapse=" ")))
	ttt=list(initial_U,initial_V,U,V,objs,objs.parts,rmse,X1)
	names(ttt)=c("initial_U","initial_V","U","V","objs","objs.parts","rmse","X1")
	return(ttt)
}

#################

#################add fixed rows of P,basically the same as qnmf_indisS_all_revise_addP_old
qnmf_indisS_all_revise_addP<-function(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,theta,qq, iter, epslog,mscale, cnormalize, addP, P_fix, tProp){
##########This is for qua matrix factorization, indiS is an indicator matrix equal to 1-NMF_indi_all; indiS_method has two options, either project descent, or gradient descent; theta is the penalty on similarity of U and NMF_indi_all matrix; qq is the value controling the updating scheme; 
##X1: original matrix 
##initial_U
##initial_V
##NMF_indi_all
##indiS_method
##UM
##DU
##VM
##DV
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
#library(pracma)
	K_common=NULL
	K=ncol(NMF_indi_all)
        if(mscale==1){
                X1=t(apply(X1, MARGIN=1, function(x)x/median(x)))
        }
	indiS=1-NMF_indi_all
	indiS[,which(P_fix==0)]=0*indiS[,which(P_fix==0)]
	if(is.matrix(initial_U) & is.matrix(initial_V)){##both initial U and V are given
	}else{
		if(is.matrix(initial_U)){ ##only initial U is given
			initial_V=compute.V(X1, initial_U)
		}else if(is.matrix(initial_V)){ ##only initial V is given
			initial_U=(compute.V(t(X1),initial_V))
		}else{ ##none of initial U or V is given
			MMM=1###MMM=2 is the worst, could pick between 1 and 2
			if(MMM==1){ ##gives better fitting
				NMF_indi_all1=NMF_indi_all
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),ncol(NMF_indi_all),"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=inis.list.qnmf[["initial_V"]]
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}else if(MMM==2){ ##gives better similarity to NMF_indi_all
				NMF_indi_all1=NMF_indi_all[,which(P_fix==1)]
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),K_common,"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=cbind(compute.V(X1, initial_U),t(addP[which(P_fix==0),]))
				initial_U=(compute.V(t(X1),initial_V))
			}else{
				initial_U=(compute.V(t(X1),t(addP[which(P_fix==0),])))
				initial_U=cbind(NMF_indi_all[,which(P_fix==1)],initial_U)
				initial_V=compute.V(X1, initial_U)
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}
		}
	}
	colnames(initial_U)=colnames(NMF_indi_all)
	colnames(initial_V)=colnames(NMF_indi_all)	#chang
	print(dim(initial_V))
	GL.list.qnmf = initialize_GL(X1,K,"qnmf")
	if(is.matrix(UM)){
		DU=diag(rowSums(UM))
	}else{
		UM=GL.list.qnmf[["UM"]]
		DU=diag(rowSums(UM))
	}
	if(is.matrix(VM)){
		DV=diag(rowSums(VM))
	}else{
		VM=GL.list.qnmf[["VM"]]
		DV=diag(rowSums(VM))
	}

	U=initial_U; V=initial_V;
####Check the initial residuals
	Xpred=U%*%t(V)
	print("This is the first iteration")
	print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	print("This is the second iteration")
####Check the initial residuals
	objs=NULL; objs.parts=NULL
	flag=0
	rrr=0
	eps=10^(-epslog) 
	extraOut=NULL
	extraOut_rnames=NULL

#######run for the first time
	s1=sum((X1-U%*%t(V))^2)
	s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
	s2=s21*alpha+s22*beta
	s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
	s3=s31*gamma+s32*roh
	s4=2*theta*sum(diag(t(U)%*%indiS))
	obj_tmp=s1+s2+s3+s4

	objs=c(objs,obj_tmp)
	extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
	extraOut_rnames=c(extraOut_rnames, 0)
########run for the first time
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_U = U;old_V = V;
		if(gamma>0){
			updateU="multiplicative"
		}else{
			updateU="HALS-col"
		}
		if(updateU=="multiplicative"){
			#optimize U by fixing V;
			a_U = (X1%*%V + alpha*initial_U + gamma*UM%*%U);
			b_U = (U%*%t(V)%*%V+alpha*U+ gamma*DU%*%U);
	                if(indiS_method=="prdescent"){
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                        update_U = update_U*NMF_indi_all
	                }else{
	                        b_U=b_U+theta*indiS
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                }
			U[U<eps]=eps #This value is very important! Because in multiplicative update, if you get a zero, it won't be updated at any further step;
			#update_U[update_U<eps]=eps
			U  = U * update_U;
			#Normalization
		}
		if(updateU=="HALS-col"){
		#	U[U<eps]=eps
			for(i in 1:ncol(U)){
				tmp.mat=NULL
				for(j in 1:ncol(V)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*U[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
		                if(indiS_method=="prdescent"){
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in NMF_indi_all matrix to be zero
		                }else{
                			a_U=a_U-theta*indiS[,i]
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                }
				update_U = update_U*(update_U>0)
				if(cnormalize==1){ ####
				#if(cnormalize==1 & P_fix[i]==1){ ####
					update_U[update_U<eps]=eps ####
					V[,i]=V[,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
					update_U=update_U/sqrt(sum(update_U^2)) ####
				} ###
				#update_U[update_U<eps]=eps
				U[,i]=update_U
			}
		}
		if(roh>0){
			updateV="multiplicative"
		}else{
			updateV="HALS-col"
		}
		if(updateV=="multiplicative"){
			#optimize V by fixing U;
			a_V = (t(X1)%*%U +beta*initial_V + roh*VM%*%V);
			b_V = (V%*%t(U)%*%U +beta*V+ roh*DV%*%V);
			b_V[b_V<eps]=eps
			update_V = (a_V/b_V)^qq;
			#update_V[update_V<eps]=eps
			#update_V = (a_V/b_V);
			V[V<eps]=eps
			V  = V * update_V;
			#Normalization
		}
		if(updateV=="HALS-col"){
		#	V[V<eps]=eps
			for(i in 1:ncol(V)){
				#if(P_fix[i]==1){
					tmp.mat=NULL
					for(j in 1:ncol(U)){
						if(j!=i){
							tmp=sum(U[,j]*U[,i])*V[,j]
							tmp.mat=cbind(tmp.mat, tmp)
						}
					}
					a_V=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_V[,i,drop=FALSE]
					b_V=(sum(U[,i]*U[,i])+beta)
					b_V[b_V<eps]=eps
					update_V=a_V/b_V
					update_V = update_V*(update_V>0)
	#			print(update_V)
					update_V[update_V<eps]=eps
					V[,i]=update_V
				#}
			}
		}
		
		s1=sum((X1-U%*%t(V))^2)
		s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
		s2=s21*alpha+s22*beta
		s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
		s3=s31*gamma+s32*roh
		s4=2*theta*sum(diag(t(U)%*%indiS))
		obj_tmp=s1+s2+s3+s4

		objs=c(objs,obj_tmp)
		#objs.parts=rbind(objs.parts,s1, s21,s22,s31,s32,s4)
		objs.parts=rbind(objs.parts, matrix(c(s1,s21,s22,s31,s32,s4),nrow=1))
		#print(res)
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(abs(objs[length(objs)]-objs[length(objs)-1])<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		
		Xpred=U%*%t(V)
#		print(paste("The fitting residual is", sum((Xpred-X1)^2)))
		if(rrr%%10==0){
			extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
			extraOut_rnames=c(extraOut_rnames, rrr)
			print(c(apply(cor(t(tProp),V),1,max), s1,s4  ))
		}
	}
####run for the last time
	extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
	extraOut_rnames=c(extraOut_rnames, rrr)
####run for the last time
	rownames(extraOut)=extraOut_rnames
	rmse=sum((Xpred-X1)^2)
	print(paste("This is the",rrr,"th iteraction"))
	print(paste("The objective function is", obj_tmp))
	print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are:",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
#	print(paste("The residual is ",res,collapse=" "))
	#print(paste("The residual is ",paste(c(res1,res2),collapse=" ")))
	#print(paste("The parts objectives are ",paste(zzz[-1],collapse=" ")))
	ttt=list(initial_U,initial_V,U,V,objs,objs.parts,rmse,X1,extraOut)
	names(ttt)=c("initial_U","initial_V","U","V","objs","objs.parts","rmse","X1","extraOut")
	return(ttt)
}

#################

##############cross-validation
#################add fixed rows of P,basically the same as qnmf_indisS_all_revise_addP, but theta parameter is changed.
qnmf_indisS_all_revise_addP_RR<-function(X1,initial_U,initial_V,NMF_indi_all,indiS_method,UM,VM,alpha,beta,gamma,roh,RR,qq, iter, epslog,mscale, cnormalize, addP, P_fix, tProp){
##########This is for qua matrix factorization, indiS is an indicator matrix equal to 1-NMF_indi_all; indiS_method has two options, either project descent, or gradient descent; theta is the penalty on similarity of U and NMF_indi_all matrix; qq is the value controling the updating scheme; 
##X1: original matrix 
##initial_U
##initial_V
##NMF_indi_all
##indiS_method
##UM
##DU
##VM
##DV
##alpha
##beta
##gamma
##roh
##theta
##qq
##iter
##epslog
#library(pracma)
	K_common=NULL
	K=ncol(NMF_indi_all)
        if(mscale==1){
                X1=t(apply(X1, MARGIN=1, function(x)x/median(x)))
        }
	indiS=1-NMF_indi_all
	indiS[,which(P_fix==0)]=0*indiS[,which(P_fix==0)]
	if(is.matrix(initial_U) & is.matrix(initial_V)){##both initial U and V are given
	}else{
		if(is.matrix(initial_U)){ ##only initial U is given
			initial_V=compute.V(X1, initial_U)
		}else if(is.matrix(initial_V)){ ##only initial V is given
			initial_U=(compute.V(t(X1),initial_V))
		}else{ ##none of initial U or V is given
			MMM=1###MMM=2 is the worst, could pick between 1 and 2
			if(MMM==1){ ##gives better fitting
				NMF_indi_all1=NMF_indi_all
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),ncol(NMF_indi_all),"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=inis.list.qnmf[["initial_V"]]
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}else if(MMM==2){ ##gives better similarity to NMF_indi_all
				NMF_indi_all1=NMF_indi_all[,which(P_fix==1)]
				inis.list.qnmf=initialize_mats(X1,ncol(NMF_indi_all1),K_common,"qnmf",NMF_indi_all1)
				initial_U=inis.list.qnmf[["initial_U"]]
				initial_V=cbind(compute.V(X1, initial_U),t(addP[which(P_fix==0),]))
				initial_U=(compute.V(t(X1),initial_V))
			}else{
				initial_U=(compute.V(t(X1),t(addP[which(P_fix==0),])))
				initial_U=cbind(NMF_indi_all[,which(P_fix==1)],initial_U)
				initial_V=compute.V(X1, initial_U)
				initial_V[,which(P_fix==0)]=t(addP[which(P_fix==0),])
				initial_U=(compute.V(t(X1),initial_V))
			}
		}
	}
	colnames(initial_U)=colnames(NMF_indi_all)
	colnames(initial_V)=colnames(NMF_indi_all)	#chang
	print(dim(initial_V))
	GL.list.qnmf = initialize_GL(X1,K,"qnmf")
	if(is.matrix(UM)){
		DU=diag(rowSums(UM))
	}else{
		UM=GL.list.qnmf[["UM"]]
		DU=diag(rowSums(UM))
	}
	if(is.matrix(VM)){
		DV=diag(rowSums(VM))
	}else{
		VM=GL.list.qnmf[["VM"]]
		DV=diag(rowSums(VM))
	}

	U=initial_U; V=initial_V;
####Check the initial residuals
	Xpred=U%*%t(V)
	#print("This is the first iteration")
	#print(paste("The fitting residual is", sum((Xpred-X1)^2)))
	#print("This is the second iteration")
####Check the initial residuals
	objs=NULL; objs.parts=NULL
	flag=0
	rrr=0
	eps=10^(-epslog) 
	extraOut=NULL
	extraOut_rnames=NULL

#######run for the first time
	theta=RR*sum((X1-U%*%t(V))^2)/(2*sum(diag(t(U)%*%indiS)))
	s1=sum((X1-U%*%t(V))^2)
	s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
	s2=s21*alpha+s22*beta
	s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
	s3=s31*gamma+s32*roh
	s4=2*theta*sum(diag(t(U)%*%indiS))
	obj_tmp=s1+s2+s3+s4

	objs=c(objs,obj_tmp)
	extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
	extraOut_rnames=c(extraOut_rnames, 0)
########run for the first time
	while(flag==0 & rrr<iter){
		rrr=rrr+1
		old_U = U;old_V = V;
		if(gamma>0){
			updateU="multiplicative"
		}else{
			updateU="HALS-col"
		}
		if(updateU=="multiplicative"){
			#optimize U by fixing V;
			a_U = (X1%*%V + alpha*initial_U + gamma*UM%*%U);
			b_U = (U%*%t(V)%*%V+alpha*U+ gamma*DU%*%U);
	                if(indiS_method=="prdescent"){
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                        update_U = update_U*NMF_indi_all
	                }else{
	                        b_U=b_U+theta*indiS
				b_U[b_U<eps]=eps
	                        update_U = (a_U/b_U)^qq;
	                }
			U[U<eps]=eps #This value is very important! Because in multiplicative update, if you get a zero, it won't be updated at any further step;
			#update_U[update_U<eps]=eps
			U  = U * update_U;
			#Normalization
		}
		if(updateU=="HALS-col"){
		#	U[U<eps]=eps
			for(i in 1:ncol(U)){
				tmp.mat=NULL
				for(j in 1:ncol(V)){
					if(j!=i){
						tmp=sum(V[,j]*V[,i])*U[,j]
						tmp.mat=cbind(tmp.mat, tmp)
					}
				}
				a_U=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+X1%*%V[,i,drop=FALSE]+alpha*initial_U[,i,drop=FALSE]
		                if(indiS_method=="prdescent"){
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                        update_U = update_U*(NMF_indi_all[,i]) ###force those entries that are zero in NMF_indi_all matrix to be zero
		                }else{
                			a_U=a_U-theta*indiS[,i]
					b_U=(sum(V[,i]*V[,i])+alpha)
					b_U[b_U<eps]=eps
					update_U=a_U/b_U
		                }
				update_U = update_U*(update_U>0)
				if(cnormalize==1){ ####
				#if(cnormalize==1 & P_fix[i]==1){ ####
					update_U[update_U<eps]=eps ####
					V[,i]=V[,i]*sqrt(sum(update_U^2))##############!!!!!!!!!!!???????
					update_U=update_U/sqrt(sum(update_U^2)) ####
				} ###
				#update_U[update_U<eps]=eps
				U[,i]=update_U
			}
		}
		if(roh>0){
			updateV="multiplicative"
		}else{
			updateV="HALS-col"
		}
		if(updateV=="multiplicative"){
			#optimize V by fixing U;
			a_V = (t(X1)%*%U +beta*initial_V + roh*VM%*%V);
			b_V = (V%*%t(U)%*%U +beta*V+ roh*DV%*%V);
			b_V[b_V<eps]=eps
			update_V = (a_V/b_V)^qq;
			#update_V[update_V<eps]=eps
			#update_V = (a_V/b_V);
			V[V<eps]=eps
			V  = V * update_V;
			#Normalization
		}
		if(updateV=="HALS-col"){
		#	V[V<eps]=eps
			for(i in 1:ncol(V)){
				#if(P_fix[i]==1){
					tmp.mat=NULL
					for(j in 1:ncol(U)){
						if(j!=i){
							tmp=sum(U[,j]*U[,i])*V[,j]
							tmp.mat=cbind(tmp.mat, tmp)
						}
					}
					a_V=-matrix(apply(tmp.mat,MARGIN=1, sum), ncol=1)+t(X1)%*%U[,i,drop=FALSE]+beta*initial_V[,i,drop=FALSE]
					b_V=(sum(U[,i]*U[,i])+beta)
					b_V[b_V<eps]=eps
					update_V=a_V/b_V
					update_V = update_V*(update_V>0)
	#			print(update_V)
					update_V[update_V<eps]=eps
					V[,i]=update_V
				#}
			}
		}
		
		s1=sum((X1-U%*%t(V))^2)
		s21=sum((U-initial_U)^2); s22=sum((V-initial_V)^2);
		s2=s21*alpha+s22*beta
		s31=sum(diag(t(U)%*%(DU-UM)%*%U)); s32=sum(diag(t(V)%*%(DV-VM)%*%V)); 
		s3=s31*gamma+s32*roh
		s4=2*theta*sum(diag(t(U)%*%indiS))
		obj_tmp=s1+s2+s3+s4

		objs=c(objs,obj_tmp)
		#objs.parts=rbind(objs.parts,s1, s21,s22,s31,s32,s4)
		objs.parts=rbind(objs.parts, matrix(c(s1,s21,s22,s31,s32,s4),nrow=1))
		#print(res)
		if(is.finite(sum(objs))){
			if(rrr>5){
				if(abs(objs[length(objs)]-objs[length(objs)-1])<0.01){
					flag=1
				}
			}
		}else{
			flag=1
			print("The program is killed because of singular updates")
		}
		
		Xpred=U%*%t(V)
#		print(paste("The fitting residual is", sum((Xpred-X1)^2)))
		if(rrr%%10==0){
			extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
			extraOut_rnames=c(extraOut_rnames, rrr)
	#		print(c(apply(cor(t(tProp),V),1,max), s1,s4  ))
		}
	}
####run for the last time
	extraOut=rbind(extraOut, c(apply(cor(t(tProp),V),1,max), s1,s4  ))
	extraOut_rnames=c(extraOut_rnames, rrr)
####run for the last time
	rownames(extraOut)=extraOut_rnames
	rmse=sum((Xpred-X1)^2)
	#print(paste("This is the",rrr,"th iteraction"))
	#print(paste("The objective function is", obj_tmp))
	#print(paste("The parameters alpha, beta, gamma, roh, theta,qq, iter, are:",paste(c(alpha, beta, gamma, roh, theta,qq, iter),collapse=" ")))
#	print(paste("The residual is ",res,collapse=" "))
	#print(paste("The residual is ",paste(c(res1,res2),collapse=" ")))
	#print(paste("The parts objectives are ",paste(zzz[-1],collapse=" ")))
	ttt=list(initial_U,initial_V,U,V,objs,objs.parts,rmse,X1,extraOut)
	names(ttt)=c("initial_U","initial_V","U","V","objs","objs.parts","rmse","X1","extraOut")
	return(ttt)
}

#################



##############
#'@importFrom nnls nnls
cv_tnmf<-function(X,indiS_method,NMF_indi_all,alpha,beta,gamma,roh, theta,qq,iter,epslog, nPerm, mscale,cnormalize){#### This still needs to be revised as how to cross-validate tnmf????
##X
##NMF_indi_all
##theta
##eta
##qq
##iter
##nPerm
        if(mscale==1){
                X=t(apply(X, MARGIN=1, function(x)x/median(x)))
        }
##########
	K=ncol(NMF_indi_all)
	indiS=1-NMF_indi_all
#########
	devs=NULL
        for(rrr in 1:nPerm){
                test_ID = sample(1:ncol(X),0.8*ncol(X))
                X_training=X[,test_ID]
		FM=GM=NULL
		initial_F=NULL
		initial_S=NULL
		initial_G=NULL
                ttt=tnmf_indisS_all_revise(X_training,initial_F,initial_S, initial_G,NMF_indi_all,indiS_method,FM,GM,alpha,beta,gamma,roh,theta,qq,iter,epslog,mscale)###Run the normal NMF
#     tnmf_indisS_all_revise<-function(X,initial_F,initial_S, initial_G,NMF_indi_all,indiS_method,FM,DF,GM,DG,alpha,beta,gamma,roh,theta,qq, iter, epslog){
                X_test=X[,-test_ID]
		if(is.infinite(ttt[["objs"]][length(ttt[["objs"]])]) | is.nan(ttt[["objs"]][length(ttt[["objs"]])])){
			save(X_training, X_test, file="X_nan.RData")
		}else{
	                F_pred = ttt$F
	                G_pred=NULL
	                for(j in 1:ncol(X_test)){
	                        fit = nnls::nnls(F_pred, X_test[,j])
	                        coffs=coef(fit)
	                        G_pred=rbind(G_pred, coffs)
	                }
	                devs = rbind(devs, c(sd(X_test),quantile(X_test),sum((F_pred%*%t(G_pred)-X_test)^2)))
		}
        }
	colnames(devs)=c("sd","0%","25%", "50%", "75%", "100%","diffs")
        return(devs)
}
