#' Convert Time Series Matrix to Connectivity Matrix for PH calculation
#' @description
#' \code{TStoConn} converts #time-points(T)*#ROIs(N) timeseries matrix or data frame to N*N symmetric connectivity matrix
#' for PH calculation while while removing ROIs(or voxels) whose values are all zero or all the same.
#' @param ts #time-points(T)*#ROIs(N) timeseries matrix
#' @param is.abs logical. TRUE returns a connectivity matrix with absolute correlation coefficient values. Default is FALSE.
#'
#' @returns K*K symmetric matrix (K is the number of survived ROIs (or voxels) after screening)
#'
#' @examples
#' TS <- matrix(runif(100*5,-1,1),nrow=100,ncol=5)
#' TStoConn(TS, is.abs=TRUE)
#' @export



TStoConn <- function(ts, is.abs=FALSE){
  #######################################
  ## when all ts values are 0 or have the exactly same ts's each other

  csa = colSums(abs(ts))
  mcsa = min(csa)
  idx = csa > 0 ### index which rois have zero ts
  if (sum(csa==mcsa)>1) idx = (csa > 0) & (csa > mcsa) ### index if some rois have same ts
  ts=ts[,idx]

  CorMat = stats::cor(ts)
  diag(CorMat) = 1
  CorMat = round(CorMat,10)
  if(is.abs) CorMat = abs(CorMat)

  if(sum(is.na(CorMat))!=0) print("NA exists")
  return(CorMat)
}

