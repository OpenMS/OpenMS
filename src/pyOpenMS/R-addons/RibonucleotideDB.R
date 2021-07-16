
# C++ signature: libcpp_pair[const Ribonucleotide *,const Ribonucleotide *] getRibonucleotideAlternatives(const libcpp_string & code)
RibonucleotideDB$$set("public","getRibonucleotideAlternatives",function(code){
  if(!is_scalar_character(code)) { stop("arg code wrong type") }
  res <- private$$py_obj$$getRibonucleotideAlternatives(py_builtin$$bytes(code,'utf-8'))
  res[[1]] <- Ribonucleotide$$new(res[[1]])
  res[[2]] <- Ribonucleotide$$new(res[[2]])
  return(res)
})