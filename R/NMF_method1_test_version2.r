NMF_method1_test_version2<-function(tg_list, data_ng, data_normalized, max_ES_cut=0.3,tProp,NMF_RR=10)
{
        tg_selected_R4_RR<-tg_list
        data.matrix0_s<-data_ng
        data_23_s<-data_normalized
        Rbase_selected_R4_RR<-Compute_Rbase_SVD(data.matrix0_s,tg_selected_R4_RR)
        stat_selected_R4_RR<-compute_IM_stat(tg_list_c=tg_selected_R4_RR)
        stat_selected_R4_RR_max<-apply(stat_selected_R4_RR,1,max)
        tg_I_id<-which(stat_selected_R4_RR_max==max(stat_selected_R4_RR_max))[1]
        tg_selected_id<-tg_I_id
        X_I<-data.matrix0_s[tg_selected_R4_RR[[tg_I_id]],]
        N<-0
        NMF_list_c<-list()
        N<-N+1
        NMF_list_c[[N]]<-tg_selected_R4_RR[[tg_I_id]]
        names(NMF_list_c)<- names(tg_I_id)
        NMF_self_c1<-Building_NMF_input_no_cancer(tg_data=X_I,tg_list=NMF_list_c,tg_list_add=list())
        V_I<-t(as.matrix(Rbase_selected_R4_RR[tg_I_id,]))
        V_I_n<-normalize_data2(V_I)
        V_c<-V_I_n
        ccc<-compute_CompRowspace_NN(data_23_s,module_RB=V_c,ROUNDS=min(nrow(V_I),3))
        ES_I<-RMSE_row(ccc)/RMSE_row(data_23_s)
        ES_base_c<-Compute_base_ES_score(ES_I,tg_selected_R4_RR)
        ES_max<-max(ES_base_c)
        tg_id_del<-c()
        track_N_cor<-list()
        track_N_cor[[1]]<-apply(cor(t(tProp),t(V_c)),1,max)
        tg_id<-setdiff(1:length(tg_selected_R4_RR),tg_I_id)
        while((ES_max>max_ES_cut)&(length(tg_id)>0))
        {
                N<-N+1
                tg_id<-which(abs(ES_base_c-max(ES_base_c))<0.01)
                tg_id<-setdiff(tg_id,tg_selected_id)
                if(length(tg_id)==0)
                {
                        break
                }
                tg_id_add<-tg_id[which(stat_selected_R4_RR_max[tg_id]==max(stat_selected_R4_RR_max[tg_id]))[1]]
                tg_selected_id<-c(tg_selected_id,tg_id_add)
                NMF_list_c<-select_R_base(tg_selected_R4_RR,tg_selected_id)
                names(NMF_list_c)<-names(tg_selected_R4_RR)[tg_selected_id]     #?????
                NMF_self_c1<-Building_NMF_input_no_cancer(tg_data=data.matrix0_s,tg_list=NMF_list_c,tg_list_add=list())
                qnmf_result_c_temp <- run_NMF(NMF_self_c1,RR0=NMF_RR,maxIter=20000, tProp)
                U_temp = qnmf_result_c_temp[["U"]]
                V_temp = qnmf_result_c_temp[["V"]]
                ########check correlation between: V and true Proportion
                V_c<-normalize_data2(t(V_temp))
                ccc<-compute_CompRowspace_NN(data_23_s,module_RB=V_c,ROUNDS=min(nrow(V_I),3))
                ES_I<-RMSE_row(ccc)/RMSE_row(data_23_s)
                ES_base_c_temp<-Compute_base_ES_score(ES_I,tg_selected_R4_RR)
                if(ES_base_c_temp[tg_id_add]>max_ES_cut)
                {
                        tg_selected_id<-setdiff(tg_selected_id,tg_id_add)
                        tg_id_del<-c(tg_id_del,tg_id_add)
                        if(length(tg_id_del)>0)
                        {
                                ES_base_c[tg_id_del]<-0
                        }
                        N<-N-1
                }
                else
                {
                        qnmf_result_c<-qnmf_result_c_temp
                        track_N_cor[[N]]<-qnmf_result_c[[4]]
                        V<- qnmf_result_c[["V"]]
                        ES_base_c<-ES_base_c_temp
                        if(length(tg_id_del)>0)
                        {
                                ES_base_c[tg_id_del]<-0
                        }
                        ES_max<-max(ES_base_c)
                        ictd.ccc= t(cor(t(tProp),V) ) 
                        o6_predict_ture = apply(ictd.ccc, 1, max)
                        o6_predict_ture2 = apply(ictd.ccc, 2, max)
                        print(N)
                        print(o6_predict_ture2)
                        print(apply(cor(t(tProp),t(Rbase_selected_R4_RR[tg_selected_id,])),1,max))
                }
        }
        return(list(qnmf_result_c,track_N_cor,NMF_self_c1))
}


