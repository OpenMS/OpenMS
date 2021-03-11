
SpectrumAlignmentScore$$set("public","score_1",function(spec1){
  if(!(is.R6(spec1) && class(spec1)[1]=="MSSpectrum")) stop("arg spec1 wrong type")
  private$$py_obj(r_to_py(spec1))
})

SpectrumAlignmentScore$$set("public","score_2",function(spec1,spec2){
  if(!(is.R6(spec1) && class(spec1)[1]=="MSSpectrum")) stop("arg spec1 wrong type")
  if(!(is.R6(spec2) && class(spec2)[1]=="MSSpectrum")) stop("arg spec2 wrong type")
  private$$py_obj(r_to_py(spec1),r_to_py(spec2))
})