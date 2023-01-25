brain_MST_segregation = function(connectivity, network){
  net = unique(network)
  within = NULL
  between = NULL
  
  for(n in net){
    idx = which(network==n)
    conn = connectivity[idx,idx]
    conn.b = apply(connectivity[idx,-idx],2,max)
    within = c(within,hclust(as.dist(1-conn),method="single")$height)
    between = c(between,conn.b)
  }
  W = mean(within)
  B = mean(between)
  bms = (W-B)/W
  return(bms)
}
