phantom.report.optk<-function(res,
                              phantom.summary.multk,
                              report_list,
                              path,
                              zipfile){
  if(!file.exists(path)){
    dir.create(file.path(path))
  }
  old_wd = getwd()
  setwd(path)
  zipfile_list = c() ## list of output file names to be zipped
  ## run_info file
  # run_info_file = paste(path,'/Run.info.csv',sep = '')
  run_info_file = 'Run.info.csv'
  zipfile_list = c(zipfile_list,run_info_file)
  run_time_stamp = res$phantom.config$time_stamp
  geneset_type = res$phantom.config$geneset_type
  max_number_clusters = res$phantom.config$maxncluster
  nsample = res$phantom.config$nsample
  reporting_pvalue = phantom.summary.multk$filter$PVAL
  reporting_nmin = phantom.summary.multk$filter$NMIN
  reporting_p_min_percentage = phantom.summary.multk$filter$p_min_percentage
  reporting_d_min_quantile = phantom.summary.multk$filter$d_min_quantile

  run_info = list(run_time_stamp = run_time_stamp,
                  geneset_type  =geneset_type,
                  max_number_clusters = max_number_clusters,
                  nsample = nsample,
                  reporting_pvalue = reporting_pvalue,
                  reporting_nmin = reporting_nmin,
                  reporting_p_min_percentage = reporting_p_min_percentage,
                  reporting_d_min_quantile = reporting_d_min_quantile)

  write.table(unlist(run_info),file = run_info_file,quote = FALSE,col.names = FALSE,sep = ',')

  report_geneset_id = unlist(lapply(report_list,function(x) x$geneset_id))
  report_geneset_name = unlist(lapply(report_list,function(x) x$geneset_name))
  report_geneset_optk = unlist(lapply(report_list,function(x) x$k))

  n = length(report_geneset_id)
  geneset_splited_name_list = list()
  geneset_splited_member_list = list()
  geneset_splited_maxsize = 0

  for (i in 1:n){
    id = report_geneset_id[i]
    obj = res$phantom.result[[id]]
    # cat('report:',id,'\n')
    ## output figs to pdf
    pdf_file_name = paste(report_geneset_name[i],'.pdf',sep = '')
    pdf(pdf_file_name,width = 10,height = 8,onefile = FALSE)
    plot.phantom(obj)
    dev.off()
    zipfile_list = c(zipfile_list,pdf_file_name)

    ## split gene sets and output to csv
    current.optk = obj$ncluster
    id.sig = which(phantom.summary.multk$summary$is_valid_cluster_all[[id]][[current.optk-1]]==TRUE &
                     phantom.summary.multk$summary$pval_all[[id]][[current.optk-1]]<=phantom.summary.multk$filter$PVAL)
    geneset_splited_name_list[[i]]=rep(0,length(id.sig)+1)
    geneset_splited_member_list[[i]] = list()
    if(current.optk==2){
      geneset_splited_name_list[[i]] = paste(obj$geneset_name,c(1,2),sep = '_')
      geneset_splited_member_list[[i]][[1]] = obj$geneset_stat$x.gene_name[which(obj$geneset_stat$x.cluster.group==1)]
      geneset_splited_member_list[[i]][[2]] = obj$geneset_stat$x.gene_name[which(obj$geneset_stat$x.cluster.group==2)]
      geneset_splited_maxsize = max(geneset_splited_maxsize,length(geneset_splited_member_list[[i]][[1]]),length(geneset_splited_member_list[[i]][[2]]))
    }else{
      tmp.members = c()
      for(j in 1:length(id.sig)){
        geneset_splited_name_list[[i]][j] = paste(obj$geneset_name,j,sep = '_')
        tmp.members = c(tmp.members,obj$geneset_stat$x.gene_name[which(obj$geneset_stat$x.cluster.group==id.sig[j])])
        geneset_splited_member_list[[i]][[j]] = obj$geneset_stat$x.gene_name[which(obj$geneset_stat$x.cluster.group==id.sig[j])]
        geneset_splited_maxsize = max(geneset_splited_maxsize,length(geneset_splited_member_list[[i]][[j]]))
      }
      if(length(id.sig)<current.optk){
        geneset_splited_name_list[[i]][length(id.sig)+1] = paste(obj$geneset_name,length(id.sig)+1,sep  = '_')
        geneset_splited_member_list[[i]][[length(id.sig)+1]] = setdiff(obj$geneset_stat$x.gene_name,tmp.members)
        geneset_splited_maxsize = max(geneset_splited_maxsize,length(geneset_splited_member_list[[i]][[length(id.sig)+1]]))
      }
    }
  }
  # cat(geneset_splited_maxsize,'\n')
  geneset_splited_name_vec = unlist(geneset_splited_name_list)
  geneset_splited_desp_vec = rep('Phantom',length(geneset_splited_name_vec))
  geneset_splited_member_mat = matrix('',length(geneset_splited_name_vec),geneset_splited_maxsize)
  tmp.id = 1
  for(i in 1:length(geneset_splited_member_list)){
    for(j in 1:length(geneset_splited_member_list[[i]])){
      geneset_splited_member_mat[tmp.id,1:length(geneset_splited_member_list[[i]][[j]])] = geneset_splited_member_list[[i]][[j]]
      tmp.id = tmp.id + 1
    }
  }
  geneset_splited_output = list(geneset_splited_name_vec,geneset_splited_desp_vec,geneset_splited_member_mat)
  geneset_splited_file = 'het.genesets.split.csv'
  write.table(geneset_splited_output,file = geneset_splited_file,quote = FALSE,col.names = FALSE,sep = ',',row.names = FALSE)
  # zipfile_list = c(zipfile_list,geneset_splited_file)
  # ## zip output results
  # zip(zipfile=zipfile, files=zipfile_list)
  setwd(old_wd)
}
