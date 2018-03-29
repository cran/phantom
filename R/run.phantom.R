#' Run phantom analysis
#'
#'This function allows user to run individual geneset heterogeneity analysis with phantom
#' @param data User provided time-course data loaded by load.data()
#' @param geneset_list User provided genesets list loaded by load.geneset(). Phanotm package provids four geneset lists from different resources: kegg, reactome, emory geneset and baylor modules. These genesets can be obtained with data(), e.g. data(kegg.geneset)
#' @param query_geneset The name of a geneset user wants to analysis. This geneset should be from the geneset_list designated by geneset_list parameter
#' @param ncluster The number of clusters within a geneset user wants to use to identify the heterogeneity of this geneset.
#' @param nsample The times of random sampling that is used to build the NULL distribution for parato front analysis.
#' @keywords run phantom
#' @export
#' @examples
#' ## load in the demo data in phantom package
#' data("time.course.data")
#'
#' ## store the analysis result in an object
#' \dontrun{
#' obj = run.phantom(data = time.course.data, geneset_list = reactome.geneset,
#'                   query_geneset ='REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES',
#'                   ncluster = 2, nsample = 1000)
#'}

run.phantom<-function(data=NULL,
                      geneset_list=NULL,
                      query_geneset=NULL,
                      ncluster = 2,
                      nsample = 1000){

  if(is.null(data)|is.null(geneset_list)|is.null(query_geneset)|is.null(ncluster)){
    stop('Error: one or more of the input arguments:\n(data,geneset_list,query_geneset,ncluster)\nneed to be specified!\n')
  }

  id = which(unlist(lapply(geneset_list,
                           function(x) x$geneset_name==query_geneset))==TRUE)
  if(length(id)==0){
    stop('Error: the query gene set does not exist!\n')
  }else if(length(id)>1){
    stop('Error: the query gene set matches multiple entries in the gene set file!\n')
  }

  obj = phantom(data,
                geneset_list[[id]],
                nsample = nsample,
                # run_gapstat = 0,
                # plot = 0,
                ncluster = ncluster,
                random_clusters = list())$obj
  plot.phantom(obj)
  return(obj)
}
