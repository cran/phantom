#' Load geneset function
#'
#'This function allows user to load in the genesets in a format that can be used by
#'run.phantom and run.phantom.batch functions.
#' @param geneset_file the genesets data file path and name on user's local machine
#' @keywords load geneset
#' @export

load.geneset<-function(geneset_file){
  geneset_raw_table = read.csv(geneset_file,header = FALSE)
  geneset_list_name = as.character(geneset_raw_table[,1])
  n = length(geneset_list_name)
  geneset_matrix = as.matrix(geneset_raw_table)
  geneset_list = c()
  for (i in 1:length(geneset_list_name)){
    tmp = setdiff(as.character(geneset_matrix[i,3:dim(geneset_matrix)[2]]),c(""))
    geneset_list[[i]] = list(geneset_type = geneset_file,
                             geneset_name = geneset_list_name[i],
                             geneset_member = tmp)
  }
  return(geneset_list)

}
