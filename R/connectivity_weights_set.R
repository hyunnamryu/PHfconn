#' Persistent homology-based weight sets of functional connectivity
#' @description \code{connectivity_weights_set} calculate dimension 0 (connected components) and dimension 1 (cycles) persistent homology and generate death weights of dimension 0 topological objects and birth weights of dimension 1 topological objects.
#' @param connectivity a numeric N*N matrix or data frame, of which element represents weighted connectivity between nodes. N is the number of nodes. (e.g. Pearson correlation matrix with diagonal element of 1).
#'
#' @return W : a set of all weights in connectivity matrix.
#' @return WB0: a set of added weights when the number of dimension 0 (connected components) topological features decreases by 1.
#' @return WB1: a set of added weights when the number of dimension 1 (cycles) topological features increases by 1.
#' @return WB10: a set of weights making cycles until all nodes are connected. It is a subset of WB1 and complement of WB11. i.e. WB10 U WB11 = WB1.
#' @return WB11: a set of weights making cycles after all nodes are connected. It is a subset of WB1 and complement of WB10. i.e. WB10 U WB11 = WB1.
#'
#' @examples
#' a <- runif(45,-0.5,0.8)
#' A <- vectomat(a)
#' connectivity_weights_set(A)
#' @export

connectivity_weights_set = function(connectivity){
  # W: all weights
  # WB0: weights in B0
  # WB1: weights making cycles
  # WB10: weights making cycles until all nodes are connected
  # WB11: weights making cycles after all ROIs are connected

  W = sort(round(uppermat(connectivity),10),decreasing = T)
  WB0 = round(1-hclust(as.dist(1-connectivity),method="single")$height,10)
  WB1 = W[!(W %in% WB0)]
  WB10 = WB1[WB1>=min(WB0)]
  WB11 = WB1[WB1<min(WB0)]
  out = list(W=W,WB0=WB0,WB1=WB1,WB10=WB10,WB11=WB11)
  return(out)
}


