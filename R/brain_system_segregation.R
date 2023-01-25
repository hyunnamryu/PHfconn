brain_system_segregation = function(connectivity, network){
  net = unique(network)
  within = NULL
  between = NULL
  
  for(n in net){
    idx = which(network==n)
    within = c(within,uppermat(connectivity[idx,idx]))
    between = c(between,as.numeric(connectivity[idx,-idx]))
  }
  W = mean(atanh(within))
  B = mean(atanh(between))
  bss = (W-B)/W
  return(bss)
}

