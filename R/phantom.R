#' @importFrom stats cutree dist hclust quantile
#' @import MASS
#' @import NMF
#' @import RColorBrewer
#' @import gplots
#' @useDynLib phantom, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices colorRampPalette dev.off pdf

phantom<-function(data,
                  geneset,
                  ncluster = 2,
                  nsample = 1000,
                  random_clusters){
  if(ncluster<2){
    stop('ncluster should be equal or larger than 2\n')
  }

  ## Predefine related parameters or functions
  MYCLUSTERCOLOR = c('red','seagreen','orange','blue','purple','cyan','yellow','violet','gray66','green')

  cat('Gene set: ',geneset$geneset_name,'\n',sep='')

  gene_name = rownames(data)
  data.nrow = dim(data)[1]
  data.ncol = dim(data)[2]

  match1 = setdiff(match(geneset$geneset_member,gene_name),NA)
  x = matrix(data[match1,],ncol=data.ncol)
  rownames(x) = gene_name[match1]
  colnames(x) = colnames(data)
  if(dim(x)[1]>ncluster){
    x.nrow = dim(x)[1]
    x.ncol = dim(x)[2]
    x.gene_name = gene_name[match1]


    ## cluster analysis
    x.cluster.d = dist(x, method = "euclidean")
    x.cluster.fit <- hclust(x.cluster.d, method="complete")
    x.cluster.group <- cutree(x.cluster.fit, k=ncluster)
    x.cluster.color = rep('',dim(x)[1])
    x.cluster.size = rep(NA,ncluster)
    for(i in 1:ncluster){
      x.cluster.color[which(x.cluster.group==i)] = MYCLUSTERCOLOR[i]
      x.cluster.size[i] = length(which(x.cluster.group==i))
    }

    # Calculate cluster heteogeneity parameters:d and p
    x.mean = colMeans(x)

    x.cluster.mean = matrix(0,nrow = ncluster,ncol = x.ncol)
    x.cluster.mean1 = matrix(0,nrow = ncluster,ncol = x.ncol)
    x.cluster.d = rep(0,ncluster)
    x.cluster.p = rep(0,ncluster)

    for(i in 1:ncluster){
      x.cluster.mean[i,] = colMeans(matrix(x[which(x.cluster.group==i),],ncol = x.ncol))
      x.cluster.mean1[i,] = colMeans(matrix(x[which(x.cluster.group!=i),],ncol = x.ncol))
      x.cluster.d[i] = het.cluster.dist(rbind(x.cluster.mean[i,],x.cluster.mean1[i,]))
      x.cluster.p[i] = min(length(which(x.cluster.group==i))/x.nrow,1-length(which(x.cluster.group==i))/x.nrow)
    }


    geneset_stat = list(x = x,
                        ncluster = ncluster,
                        x.gene_name = x.gene_name,
                        x.cluster.group = x.cluster.group,
                        x.cluster.color = x.cluster.color,
                        x.cluster.size = x.cluster.size,
                        x.mean = x.mean,
                        x.custer.mean = x.cluster.mean,
                        x.cluster.d = x.cluster.d,
                        x.cluster.p = x.cluster.p,
                        MYCLUSTERCOLOR = MYCLUSTERCOLOR)

    ## random sampling
    time1 = proc.time()
    obj.random = clustering.random.genesets(data = data,
                                                     nsample=nsample,
                                                     nsize=x.nrow,
                                                     ncluster=ncluster,
                                                     random_clusters = random_clusters)
    random_geneset_stat = obj.random$obj
    random_clusters = obj.random$random_clusters
    time2 = proc.time()
    time.elapse = time2-time1
    cat('Time used for random sampling:',time.elapse[3],'\n')
    ## Pareto front analysis
    time1 = proc.time()

    pf = paretoFrontTest2C(x.cluster.d,
                            x.cluster.p,
                            c(1:ncluster),
                            x.cluster.size,
                            random_geneset_stat$x.cluster.d.null.list,
                            random_geneset_stat$x.cluster.p.null.list,
                            random_geneset_stat$x.cluster.rsid.null.list)
    pf$color = colorRampPalette(brewer.pal(11,"Spectral"))(pf$nlevel)
    pf$cluster.color = MYCLUSTERCOLOR[pf$cluster.id]
    time2 = proc.time()
    time.elapse = time2-time1
    cat('Time used for pareto front test:',time.elapse[3],'\n')
    pval = pf$pval
    cat(pval,'\n')


    run_param = list(ncluster = ncluster,
                     nsample = nsample)
    ## save result
    obj = list(geneset_type = geneset$geneset_type,
               geneset_name = geneset$geneset_name,
               geneset_stat = geneset_stat,
               random_geneset_stat = random_geneset_stat,
               run_param = run_param,
               pf = pf,
               status = 'TEST',
               ncluster = ncluster
    )

  }else{
    ## save result
    obj = list(geneset_name = geneset$geneset_name,
               status = 'NO TEST',
               ncluster = ncluster
    )
  }

  return(list(obj = obj,
              random_clusters = random_clusters))

}
