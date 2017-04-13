phantom.choose.optk<-function(phantom.result.allk,
                              PVAL,
                              NMIN,
                              p_min_percentage = 0.05,
                              d_min_quantile = 0.05
                              ){
  ### summarize runs for multiple k (from k=2 to maxncluster)
  ngeneset = length(phantom.result.allk[[1]]$phantom.result)
  geneset_name_list = rep('',ngeneset)
  for (i in 1:ngeneset){
    geneset_name_list[i] = phantom.result.allk[[1]]$phantom.result[[i]]$geneset_name
  }
  nk = length(phantom.result.allk)
  k_all = c(1:nk)+1
  pval = matrix(NA,ngeneset,nk)
  pval_all = list()
  is_new_cluster_all = list()
  has_new_valid_sig_cluster = matrix(NA,ngeneset,nk)
  is_valid_cluster_all = list()
  is_valid_new_cluster_all = list()
  p_all = list()
  p_null_all = list()
  d_all = list()
  d_null_all = list()
  cluster_size_all = list()
  eff_size_all = list()
  cluster_color_all = list()
  cluster_pch_all = list()
  nsig_cluster = matrix(0,ngeneset,nk)
  ## temporary parameters
  p = matrix(NA,ngeneset,nk)
  csize = matrix(NA,ngeneset,nk)
  cluster.eff.size = matrix(NA,ngeneset,nk)
  
  for (i in 1:ngeneset){
    pval_all[[i]] = list()
    is_new_cluster_all[[i]] = list()
    is_valid_cluster_all[[i]] = list()
    is_valid_new_cluster_all[[i]] = list()
    p_all[[i]] = list()
    p_null_all[[i]] = list()
    d_all[[i]] = list()
    d_null_all[[i]] = list()
    cluster_size_all[[i]] = list()
    eff_size_all[[i]] = list()
    cluster_color_all[[i]] = list()
    cluster_pch_all[[i]] = list()
    for (j in 1:nk){
      if(phantom.result.allk[[j]]$phantom.result[[i]]$status=="TEST"){
        pval[i,j] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$pval
        
        pval_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$pval.obj.all[3,]
        d_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$x
        d_null_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$x.null
        p_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$y
        p_null_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$y.null
        
        cluster_size_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$pval.obj.all[4,]
        eff_size_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$pval.obj.all[5,]
        cluster_color_all[[i]][[j]] = phantom.result.allk[[j]]$phantom.result[[i]]$pf$cluster.color
        
        ## if the sub-cluster is a new cluster for current k
        if(j==1){
          is_new_cluster_all[[i]][[j]] = rep(TRUE,2)
          # cluster_pch_all[[i]][[j]] = rep(19,2)
        }else{
          group_old = phantom.result.allk[[j-1]]$phantom.result[[i]]$geneset_stat$x.cluster.group
          group_current = phantom.result.allk[[j]]$phantom.result[[i]]$geneset_stat$x.cluster.group
          K_current = length(unique(group_current))
          is_new_cluster_all[[i]][[j]] = rep(TRUE,K_current)
          # cluster_pch_all[[i]][[j]] = rep(19,K_current)
          
          for(k in 1:K_current){
            if(length(unique(group_old[which(group_current==k)]))==1 & length(unique(group_old[which(group_current!=k)]))==(K_current-2)){
              is_new_cluster_all[[i]][[j]][k] = FALSE
              # cluster_pch_all[[i]][[j]][k] = 4 ## repetitive clusters
            }
          }
        }
        
        ## if the sub-cluster is a valid cluster 
        ## (meets two criteria: p_min_percentage/NMIN, and d_min_quantile)
        is_valid_cluster_all[[i]][[j]] = (d_all[[i]][[j]]>=quantile(d_null_all[[i]][[j]],d_min_quantile) &
                                            eff_size_all[[i]][[j]]>=NMIN & 
                                            p_all[[i]][[j]]>=p_min_percentage)
        
        is_valid_new_cluster_all[[i]][[j]] = is_valid_cluster_all[[i]][[j]] & is_new_cluster_all[[i]][[j]]
        
        has_new_valid_sig_cluster[i,j] = (length(which(is_valid_new_cluster_all[[i]][[j]]==TRUE &
                                                         pval_all[[i]][[j]]<=PVAL))>0)
        
        cluster_pch_all[[i]][[j]] = rep(19,length(is_valid_cluster_all[[i]][[j]]))
        cluster_pch_all[[i]][[j]][which(is_valid_new_cluster_all[[i]][[j]]==FALSE)] = 4
        
        #       nsig_cluster[i,j] = length(which(pval_all[[i]][[j]]<=PVAL & 
        #                                          eff_size_all[[i]][[j]]>=NMIN & 
        #                                          # is_new_cluster_all[[i]][[j]]==1 &
        #                                          d_all[[i]][[j]]>=quantile(d_null_all[[i]][[j]],d_min_quantile) &
        #                                          p_all[[i]][[j]]>=quantile(p_null_all[[i]][[j]],p_min_quantile)))
        nsig_cluster[i,j] = length(which(pval_all[[i]][[j]]<=PVAL & 
                                           is_valid_cluster_all[[i]][[j]]==TRUE))
        if(j==1){
          nsig_cluster[i,j] = min(1,nsig_cluster[i,j])
        }
        
      }
    }
  }
  
  
  phantom.summary.multk = list()
  phantom.summary.multk$filter = list(PVAL = PVAL,
                                      NMIN = NMIN,
                                      p_min_percentage = p_min_percentage,
                                      d_min_quantile = d_min_quantile)
  
  phantom.summary.multk$summary = list(ngeneset = ngeneset,
                                       geneset_name_list = geneset_name_list,
                                       nk = nk,
                                       k_all = k_all,
                                       pval = pval,
                                       pval_all = pval_all,
                                       is_new_cluster_all = is_new_cluster_all,
                                       has_new_valid_sig_cluster = has_new_valid_sig_cluster,
                                       is_valid_cluster_all = is_valid_cluster_all,
                                       is_valid_new_cluster_all = is_valid_cluster_all,
                                       p_all = p_all,
                                       p_null_all = p_null_all,
                                       d_all = d_all,
                                       d_null_all = d_null_all,
                                       cluster_size_all = cluster_size_all,
                                       eff_size_all = eff_size_all,
                                       cluster_color_all = cluster_color_all,
                                       cluster_pch_all = cluster_pch_all,
                                       nsig_cluster = nsig_cluster)
  ## choose gene sets (k) to report
  report_list = list()
  report_list_n = 0
  for(id in 1:length(phantom.summary.multk$summary$geneset_name_list)){
    tmp_k_id = which(phantom.summary.multk$summary$has_new_valid_sig_cluster[id,]==TRUE)
    if(length(tmp_k_id)>0){
#       tmp_nsig_cluster = nsig_cluster[id,tmp_k_id]
#       tmp_opt_k_id = which(tmp_nsig_cluster==max(tmp_nsig_cluster))[1]
#       if(length(tmp_opt_k_id)>0){
#         report_list_n = report_list_n + 1
#         report_list[[report_list_n]] = list(geneset_id = id,
#                                             geneset_name = phantom.summary.multk$summary$geneset_name_list[id],
#                                             k = phantom.summary.multk$summary$k_all[tmp_k_id[tmp_opt_k_id]])
#       }

      report_list_n = report_list_n + 1
      report_list[[report_list_n]] = list(geneset_id = id,
                                          geneset_name = phantom.summary.multk$summary$geneset_name_list[id],
                                          k = phantom.summary.multk$summary$k_all[tmp_k_id[1]])
    }
  }
  
  ## official plots to report (optimal k, significant sub-clusters)
  report_geneset_id = unlist(lapply(report_list,function(x) x$geneset_id))
  report_geneset_optk = unlist(lapply(report_list,function(x) x$k))
  
  phantom.result.optk = list()
  for (i in 1:ngeneset){
    id = which(report_geneset_id==i)
    if(length(id)>0){
      phantom.result.optk[[i]] = phantom.result.allk[[report_geneset_optk[id]-1]]$phantom.result[[i]]
    }else{
      ## if none of the k tests return significant result, record k=2 for reporting
      phantom.result.optk[[i]] = phantom.result.allk[[1]]$phantom.result[[i]]
    }
  }
  
  obj = list(phantom.result.optk = phantom.result.optk,
             phantom.summary.multk = phantom.summary.multk,
             report_list = report_list)
  
}