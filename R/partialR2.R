partialR2 = function(reduced.lm, full.lm){
  SSE.R = sum(reduced.lm$residuals^2)
  SSE.F = sum(full.lm$residuals^2)
  pr2 = (SSE.R-SSE.F)/SSE.R
  return(pr2)
}