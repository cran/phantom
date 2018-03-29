#' Print all geneset names in geneset list
#'
#'This function allows user to print all the geneset names in the loaded geneset list such that user can find the specific query_geneset name for run.phantom funciton.
#' @param geneset_list User provided genesets list loaded by load.geneset(). Phanotm package provids four geneset lists from different resources: kegg, reactome, emory geneset and baylor modules. These genesets can be obtained with data(), e.g. data(kegg.geneset)
#' @keywords print geneset names
#' @export
#' @examples
#'
#' ## store all the geneset names in one vector
#'  \dontrun{
#' g.names = geneset.names(reactome.geneset)
#' }

geneset.names = function(geneset_list = NULL){
  geneset_names = c()
  for(i in 1:length(geneset_list)){
    geneset_names = c(geneset_names,geneset_list[[i]][[2]])
  }
  return(geneset_names)
}


