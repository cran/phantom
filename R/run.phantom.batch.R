#' Run phantom batch analysis
#'
#'This function allows user to run batch analysis of a full geneset list with phantom, and download the identified heterogeneous genesets
#' @param data User provided time-course data loaded by load.data()
#' @param geneset_list User provided genesets list loaded by load.geneset(). Phanotm package provids four geneset lists from different resources: kegg, reactome, emory geneset and baylor modules. These genesets can be obtained with data(), e.g. data(kegg.geneset)
#' @param maxncluster The maximum number of clusters within a geneset user wants to test with. All numbers from 1 to ncluster will be tested and an optimal cluster number will be selected to identify the heterogeneity of this geneset.
#' @param nsample The times of random sampling that is used to build the NULL distribution for parato front analysis.
#' @param report_pval The maximum p value of a geneset that will be reported as a significant heterogeneous geneset. Genesets with p value larger than report_pval wil not be reported
#' @param report_nmin The minmum size of subcluster in a geneset that will be reported as a significant heterogeneous geneset.
#' @param output_dir The directory where user wants to put the phantom batch analysis results
#' @keywords run phantom batch
#' @export
#' @examples
#' ## load in the demo data and geneset in phantom package
#' data("time.course.data")
#' data("kegg.geneset")
#'
#' ## store the analysis result in an object
#' \dontrun{obj = run.phantom.batch(data = time.course.data, geneset_list = kegg.geneset,
#'                   maxncluster = 5, nsample = 1000, report_pval = 0.05, report_nmin = 5,
#'                   output_dir = file.path(getwd(),'/phantom_result'))}


run.phantom.batch<-function(data = NULL,
                            geneset_list = NULL,
                            maxncluster = 5,
                            nsample = 1000,
                            report_pval = 0.05,
                            report_nmin = 5,
                            output_dir = file.path('./phantom_result')
){
  if(!file.exists(output_dir)){
    dir.create(output_dir)
  }


  ## configuration object
  time_stamp = gsub(":","-",as.character(Sys.time()),fixed = T)
  phantom_run_prefix = paste('Run_',time_stamp,sep = '')
  # phantom.pdf.name = paste('Run_',time_stamp,'.pdf',sep = '')
  phantom.zipfile.name = paste('Run_',time_stamp,'.zip',sep = '')
  phantom.config = list(#geneset_type = geneset_file,
                        maxncluster = maxncluster,
                        nsample = nsample,
                        time_stamp = time_stamp,
                        phantom.zipfile.name = phantom.zipfile.name,
                        phantom_run_prefix = phantom_run_prefix)

  ## Run Phantom
  phantom.result.allk = list()
  random_clusters = list()

  for(id in 1:length(geneset_list)){
    for(k in 2:maxncluster){
      set.seed(1)
      cat('k=',k,'id=',id,'\n')

      if(id==1){
        phantom.result.allk[[k-1]] = list(ncluster = k,
                                          phantom.result = list())
      }

      phantom.obj = phantom(data = data,
                            geneset = geneset_list[[id]],
                            ncluster = k,
                            nsample = nsample,
                            random_clusters = random_clusters)
      random_clusters = phantom.obj$random_clusters
      phantom.result.allk[[k-1]]$phantom.result[[id]] = phantom.obj$obj
    }
    random_clusters = list()
  }

  phantom.optk.obj = phantom.choose.optk(phantom.result.allk,
                                         PVAL = report_pval,
                                         NMIN = report_nmin)

  res = list(phantom.config = phantom.config,
             phantom.result = phantom.optk.obj$phantom.result.optk,
             phantom.result.allk = phantom.result.allk)

  report_file_path = file.path(output_dir,phantom.config$phantom_run_prefix)
  phantom.report.optk(res,
                      phantom.optk.obj$phantom.summary.multk,
                      phantom.optk.obj$report_list,
                      path = report_file_path,
                      zipfile = phantom.config$phantom.zipfile.name)
  return(res)
}
