#' Persistent homology-based weight sets of functional connectivity
#' @description \code{connectivity_weights_set} calculate dimension 0 (connected components) and dimension 1 (cycles) persistent homology and generate death weights of dimension 0 topological objects and birth weights of dimension 1 topological objects.
#' @param connectivity a numeric N*N matrix or data frame, of which element represents weighted connectivity between nodes. N is the number of nodes. (e.g. Pearson correlation matrix with diagonal element of 1).
#' @param na.rm logical; whether remove NA values in connectivity matrix
#'
#' @return W : a set of all weights in connectivity matrix.
#' @return WB0: a set of added weights when the number of dimension 0 (connected components) topological features decreases by 1.
#' @return WB1: a set of added weights when the number of dimension 1 (cycles) topological features increases by 1.
#' @return WB10: a set of weights making cycles until all nodes are connected. It is a subset of WB1 and complement of WB11. i.e. WB10 U WB11 = WB1.
#' @return WB11: a set of weights making cycles after all nodes are connected. It is a subset of WB1 and complement of WB10. i.e. WB10 U WB11 = WB1.
#'
#' @references
#' Ryu, H., Habeck, C., Stern, Y., & Lee, S. (2023). Persistent homology-based functional connectivity and its association with cognitive ability: Life-span study. Human Brain Mapping, 44( 9), 3669– 3683. https://doi.org/10.1002/hbm.26304
#'
#' @examples
#' a <- runif(45,-0.5,0.8)
#' A <- vectomat(a)
#' connectivity_weights_set(A)
#' @export

connectivity_weights_set = function(connectivity, na.rm=FALSE){
  # W: all weights
  # WB0: weights in B0
  # WB1: weights making cycles
  # WB10: weights making cycles until all nodes are connected
  # WB11: weights making cycles after all ROIs are connected

  if(na.rm) {
    idx = !(1:ncol(connectivity) %in% which(apply(connectivity,1,function(x) sum(!(is.na(x)))==1)))
    connectivity=connectivity[idx,idx]
  }

  #### remove totally disconnected nodes
  ## assumption: complete graph; colSums(connectivity)=0 means the node are disconnected to all the other nodes
  idx.connected = which(colSums(connectivity)!=0 & !(is.na(colSums(connectivity))))
  connectivity = connectivity[idx.connected,idx.connected]

  W = sort(round(uppermat(connectivity),10),decreasing = T)
  WB0 = round(1-stats::hclust(stats::as.dist(1-connectivity),method="single")$height,10)
  WB1 = W[!(W %in% WB0)]
  WB10 = WB1[WB1>=min(WB0)]
  WB11 = WB1[WB1<min(WB0)]
  out = list(W=W,WB0=WB0,WB1=WB1,WB10=WB10,WB11=WB11)
  return(out)
}

