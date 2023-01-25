library(igraph)
library(brainGraph)



graph_topo = function(connectivity, p=0.1){
  N = nrow(connectivity)
  M = N*(N-1)/2 ## number of possible edges
  N.E = floor(p*M) ## 10% edges
  E.thr = sort(uppermat(connectivity),decreasing = T)[N.E]
  
  binMat = connectivity>E.thr
  gg = graph_from_adjacency_matrix(binMat,mode = "undirected",diag=F)
  
  ####### measure of integration: 
  ## chracteristic path length / com
  CPL = mean_distance(gg)
  ## global efficiency
  # GE = 2*sum(1/distances(gg)[upper.tri(distances(gg))])/N/(N-1)
  GE = brainGraph::efficiency(gg,type="global")
  
  
  ######## measure of segregation:
  ## transivitity
  TV = transitivity(gg)
  ## local efficiency
  # LE = mean(brainGraph::efficiency(gg,type="local"))
  
  out = list(CPL=CPL,TV=TV,GE=GE)
  return(out)
}


