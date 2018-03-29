#' Load data function
#'
#'This function allows user to load in the time-course gene differential expression data in a format that can be used by
#'run.phantom and run.phantom.batch functions.
#' @param filename the time-course gene differential expression data file path and name on user's local machine
#' @keywords load data
#' @export

load.data<-function(filename){
  ## load in test statistic of time-course/longitudinal data
  df_tstat = read.csv(filename)
  gene_name = df_tstat$gene_symbol
  data = as.matrix(df_tstat[,2:dim(df_tstat)[2]])
  data.nrow = dim(data)[1]
  data.ncol = dim(data)[2]
  rownames(data) = gene_name

  return (data)
}
