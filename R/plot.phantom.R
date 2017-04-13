

plot.phantom<-function(obj){
  ## extract data from object
  geneset_type = obj$geneset_type
  geneset_name = obj$geneset_name
  x = obj$geneset_stat$x
  x.ncol = dim(x)[2]
  ncluster = obj$geneset_stat$ncluster
  x.cluster.group = obj$geneset_stat$x.cluster.group
  x.cluster.color = obj$geneset_stat$x.cluster.color
  x.mean = obj$geneset_stat$x.mean
  x.cluster.mean = obj$geneset_stat$x.custer.mean
  x.cluster.d = obj$geneset_stat$x.cluster.d
  x.cluster.p = obj$geneset_stat$x.cluster.p

  pf = obj$pf


    ## split geneset name if its too long
    MAXNCHAR = 60
    geneset_name_nchar = nchar(geneset_name)
    geneset_output_name = c('')
    tmp.N = floor(geneset_name_nchar/MAXNCHAR)
    if (tmp.N<1){
      geneset_output_name = paste(geneset_name,'\n')
    }else{
      for (tmp.i in 1:(tmp.N+1)){
        if (tmp.i==1){
          geneset_output_name = paste(geneset_output_name,substr(geneset_name,1+(tmp.i-1)*MAXNCHAR,tmp.i*MAXNCHAR),sep = '')
        }else if(tmp.i==(tmp.N+1)){
          geneset_output_name = paste(geneset_output_name,'\n',substr(geneset_name,1+(tmp.i-1)*MAXNCHAR,nchar(geneset_name)),'\n',sep = '')
        }else{
          geneset_output_name = paste(geneset_output_name,'\n',substr(geneset_name,1+(tmp.i-1)*MAXNCHAR,tmp.i*MAXNCHAR),sep = '')
        }
      }
    }

  # par(oma = c(1,2,3,1),mar = c(6,5,4,1))
#   layout(matrix(c(rep(1,5),rep(2,3),rep(1,4),c(3,4,4,4)),nrow = 2,byrow = TRUE))
#   hist(rnorm(100))
  ## Plot data and statistics
  par(oma = c(1,2,3,1),mar = c(6,6,5,1))
  # layout(matrix(c(1,1,2,3),2,2))
  # layout(matrix(c(rep(1,5),rep(2,3),rep(1,4),c(3,4,4,4)),nrow = 2,byrow = TRUE))
  # layout(matrix(c(rep(1,6),rep(2,4),rep(1,5),c(3,rep(4,4))),nrow = 2,byrow = TRUE))
  layout(matrix(c(rep(1,6),rep(2,4),rep(1,4),3,rep(4,5)),nrow = 2,byrow = TRUE))


  ## Plot heatmap
  MINX = min(x)
  MAXX = max(x)
  MYX = max(abs(MINX),abs(MAXX))
  max.x = max(x)
  aheatmap(obj$geneset_stat$x,
           color = bluered(10),
           Colv = NA,
           scale = 'none',
           breaks = 0,
           # labCol = NA,
           # annRow = list(Clusters = paste('sub-cluster',x.cluster.group)),
           annRow = list(Clusters = factor(paste('sub-cluster',x.cluster.group),levels = paste('sub-cluster',1:max(x.cluster.group)))),
           annColors = list(Clusters = obj$geneset_stat$MYCLUSTERCOLOR),
           main = 'Fig.1 Heatmap',
           fontsize = 10,
           sub = geneset_output_name)

#   ##
#   boxplot(x,ylim = c(MINX,MAXX))
#
#   for(i in 1:ncluster){
#     if(i==1){
#       boxplot(matrix(x[which(x.cluster.group==i),],ncol = x.ncol),
#               col = x.cluster.color[which(x.cluster.group==i)][1],
#               ylim = c(MINX,MAXX))
#     }else{
#       boxplot(matrix(x[which(x.cluster.group==i),],ncol = x.ncol),
#               col = x.cluster.color[which(x.cluster.group==i)][1],
#               ylim = c(MINX,MAXX),add = TRUE)
#     }
#
#   }
#
  ## plot average t score
  plot(colMeans(x),
       ylim = c(MINX,MAXX),
       type = 'l',
       mgp=c(2,1,0),
       xlab = 'Time points',
       ylab = 'Mean of test statistic')
  points(colMeans(x),ylim = c(MINX,MAXX),pch = 19)
  for(i in 1:ncluster){
    lines(colMeans(matrix(x[which(x.cluster.group==i),],ncol = x.ncol)),
          col = x.cluster.color[which(x.cluster.group==i)][1])
    points(colMeans(matrix(x[which(x.cluster.group==i),],ncol = x.ncol)),
           col = x.cluster.color[which(x.cluster.group==i)][1],
           pch = 19)
  }
  title('Fig.2 Change by sub-clusters', line = 1,cex.main = 1)
  # Plot Pareto front test result
  plot.pareto.front(obj)


#   title(geneset_output_name,outer = TRUE)
  title('Phantom plots',outer = TRUE,cex.main = 2)

}
