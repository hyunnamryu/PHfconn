PH_meas = function(connectivity){
  #### remove totally disconnected nodes
  ## assumption: complete graph; colSums(connectivity)=0 means the node are disconnected to all the other nodes
  idx.connected = which(colSums(connectivity)!=0)
  connectivity = connectivity[idx.connected,idx.connected]
  
  W_out = connectivity_weights_set(connectivity)
  # BC = Betti_curve(connectivity)
  
  N = dim(connectivity)[1]
  M = length(W_out$WB0) # number of edges in MST ; N-1
  L = length(W_out$WB1) # number of total cycles; N(N-1)/2-(N-1) = (N-2)(N-1)/2
  
  # mean weights of tree(MST); global weight (cost) efficiency of tree (<-> characteristic path length)
  mB0 = mean(W_out$WB0)
  
  # range of weights of tree(MST); global weight (cost) efficiency of tree (<-> diameter)
  rB0 = diff(range(W_out$WB0))
  
  # AUC under B1 curve until all nodes are connected
  mB1 = mean(W_out$WB10)
  
  # height of B1 curve until all nodes are connected - the total number of cycles when all nodes are connected
  hB1 = length(W_out$WB10)/L
  #hB1 = mean(as.numeric(lapply(W_out$WB0,function(x) sum(W_out$W>=x)))/(M+L))
  
  out = list(BS=mB0,BD=rB0,CS=mB1,CF=hB1)
  
  return(out)
}

