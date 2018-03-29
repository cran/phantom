

het.cluster.dist<-function(x,# matrix with two rows
                           method = 'euclidean'){
  if(method=='euclidean'){
    d = sqrt(sum((x[1,]-x[2,])^2))
  }
  return(d)
}
