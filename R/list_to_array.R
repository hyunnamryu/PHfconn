
list_to_array=function(lobj){
  out = array(unlist(lobj),dim=dim(lobj),dimnames=dimnames(lobj))
  return(t(out))
}


