#' Convert Vector to symmetric Matrix
#' @description \code{vectomat} converts vector with N*(N-1)/2 values to N*N symmetric matrix.
#' @param vec a vector with N*(N-1)/2 numeric values
#' @param val.diag a value which is assigned to the diagonal. Default is 1.
#'
#' @returns N*N symmetric matrix
#'
#' @examples
#' a <- runif(6,-0.5,0.8)
#' vectomat(a)
#' @export


vectomat = function(vec,val.diag=1){
  N=ceiling(sqrt(2*length(vec)))
  mat = matrix(0,ncol=N,nrow=N)
  mat[upper.tri(mat)] = vec
  mat = mat+t(mat)
  diag(mat)=val.diag
  return(mat)
}
