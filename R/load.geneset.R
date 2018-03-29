#' Load geneset in gmt file format
#'
#'This function allows user to load in the genesets that are downloaded from GSEA (gmt format)
#'into a format that can be used by run.phantom and run.phantom.batch functions.
#' @param geneset_file the genesets gmt file path and name on user's local machine
#' @keywords load geneset
#' @export

load.geneset<-function(geneset_file){
  # geneset_raw_table = read.csv(geneset_file,header = FALSE)
  geneset_raw_table = read.gmt(geneset_file)
  geneset_list_name = as.character(names(geneset_raw_table))
  n = length(geneset_list_name)
  geneset_list = c()
  for (i in 1:length(geneset_list_name)){
    geneset_list[[i]] = list(geneset_type = geneset_file,
                             geneset_name = geneset_list_name[i],
                             geneset_member = geneset_raw_table[[i]])
  }
  return(geneset_list)
}
