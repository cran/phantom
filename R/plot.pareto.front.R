plot.pareto.front<-function(obj,
                            plot.colorbar = TRUE){
  pf = obj$pf
  xmin = min(pf$x,pf$x.null)
  xmax = max(pf$x,pf$x.null)
#   ymin = min(pf$y,pf$y.null)
#   ymax = max(pf$y,pf$y.null)
  
  f1 = pf$x.null
  f2 = pf$y.null
  d = data.frame(f1,f2)
  
  ## color bar
  if(plot.colorbar==TRUE){
    lut = pf$color
    color.bar(lut, 0, 1, title='p-values', nticks=5)
  }
  
  
  ## pareto front
  plot(1,type = 'n',
       xlim = c(xmin,xmax),
       # ylim = c(0,0.5),
       ylim = c(0,0.5),
       mgp=c(2,1,0),
       # log = 'x',
       xlab = 'Distance of sub-cluster: d',
       ylab = 'Size of sub-cluster: p'
       )
  
  for(l in 1:pf$nlevel){
    front = d[which(pf$index == l),]
    front.sorted = front[order(front$f2,front$f1,decreasing=TRUE),]
    points(front$f1,front$f2,pch=19,col = pf$color[l],cex = 1.25)
  }
  
  tmp.id = c()
  for(i in 1:length(pf$x)){
    tmp.id = c(tmp.id,which((pf$x[i]>pf$x.null & pf$y[i]>=pf$y.null)
                            | (pf$x[i]>=pf$x.null & pf$y[i]>pf$y.null)))
  }
  
  if(length(tmp.id)>0){
    # pf.id =  min(pf$index[tmp.id])
    pf.id = min(which(pf$ref.pvals>=0.05))
    front = d[which(pf$index == pf.id),]
    front.sorted = front[order(front$f2,front$f1,decreasing=TRUE),]
    points(front.sorted$f1,front.sorted$f2,pch=19,col = 'black',cex = 1.5)
    lines(front.sorted$f1,front.sorted$f2,pch=19,lwd = 3,lty = 3,col = 'black')
  }
  for(i in 1:length(pf$x)){
    if(obj$geneset_stat$ncluster==2 & i==2){
      points(pf$x[i],pf$y[i],col = pf$cluster.color[i],pch = 15,cex = 2)
    }else{
      points(pf$x[i],pf$y[i],col = pf$cluster.color[i],pch = 15,cex = 3)
    }
    
  }
  # points(pf$x,pf$y,col = 'cyan',pch = 17,cex = 2)
  legend('topright',legend = paste('p-value:',as.character(round(pf$pval,3))),cex = 1.5)
  
  title('Fig.3 Pareto front plot',line = 1,cex.main = 1)
  # plot(f1,f2)

}