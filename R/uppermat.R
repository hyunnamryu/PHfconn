#' Extract upper diagonal elements from a symmetric matrix
#' @description \code{uppermat} extract a vector with N*(N-1)/2 upper diagonal values from a N*N symmetric matrix.
#' @param mat symmetric matrix
#'
#' @returns a numeric vector
#'
#' @examples
#' a <- runif(6,-0.5,0.8)
#' A <- vectomat(a)
#' uppermat(A)
#'
#' @export


uppermat = function(mat){
  uppermat = mat[upper.tri(mat)]
  return(uppermat)
}
