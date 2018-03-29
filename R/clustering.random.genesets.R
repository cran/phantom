
clustering.random.genesets<-function(data,## entire time-course data set
                                     nsample,## number of random genesets/modules to draw from entire dataset
                                     nsize,## size of genesets/module to draw
                                     ncluster,## number of clusters
                                     random_clusters
                                     ){
  x.cluster.rsid.null.list = rep(NA,ncluster*nsample)
  x.cluster.d.null.list = rep(NA,ncluster*nsample)
  x.cluster.p.null.list = rep(NA,ncluster*nsample)

  N = dim(data)[1]
  M = dim(data)[2]

  ## random sample genesets
  if(length(random_clusters)==0){
    x.null.list = list()
    x.cluster.fit.null.list = list()
    x.cluster.group.null.list = list()
    time1 = proc.time()
    for(j in 1:nsample){
      # id.null = sample(N)[1:nsize]
      id.null = sample(N,nsize)
      x.null = matrix(data[id.null,],ncol=M)

      x.cluster.dist.null = dist(x.null, method = "euclidean")
      x.cluster.fit.null <- hclust(x.cluster.dist.null, method="complete")
      x.cluster.group.null <- cutree(x.cluster.fit.null, k=ncluster)

      x.null.list[[j]] = x.null
      x.cluster.fit.null.list[[j]] = x.cluster.fit.null
      x.cluster.group.null.list[[j]] = x.cluster.group.null
    }
    time2 = proc.time()
    time.elapse = time2-time1
    cat('Time used for clustering random samples:',time.elapse[3],'\n')

    random_clusters = list(x.null.list = x.null.list,
                                       x.cluster.fit.null.list = x.cluster.fit.null.list)
  }else{
    cat('Re-using random samples for size: ',nsize,'\n')

    x.null.list = random_clusters$x.null.list
    x.cluster.fit.null.list = random_clusters$x.cluster.fit.null.list
    x.cluster.group.null.list = list()

    time1 = proc.time()
    for(j in 1:nsample){
#       id.null = sample(N,nsize)
#       x.null = matrix(data[id.null,],ncol=M)
#
#       x.cluster.dist.null = dist(x.null, method = "euclidean")
#       x.cluster.fit.null <- hclust(x.cluster.dist.null, method="complete")
      x.null = x.null.list[[j]]
      x.cluster.fit.null = x.cluster.fit.null.list[[j]]
      x.cluster.group.null <- cutree(x.cluster.fit.null, k=ncluster)
      x.cluster.group.null.list[[j]] = x.cluster.group.null
    }
    time2 = proc.time()
    time.elapse = time2-time1
    cat('Time used for clustering random samples:',time.elapse[3],'\n')
  }


  ## get heterogeneity parameters from random samples
  res = getNullHetParamC(x.null.list,x.cluster.group.null.list,ncluster)

  obj = list(x.cluster.d.null.list = unlist(res$d_list),
             x.cluster.p.null.list = unlist(res$p_list),
             x.cluster.rsid.null.list = unlist(res$rsid_list)
             )
  return(list(obj = obj,random_clusters = random_clusters))

}
