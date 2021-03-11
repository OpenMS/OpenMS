
SpectrumAlignment$$set("public","getSpectrumAlignment",function(result,spec1,spec2){
  if(!is_list(result)) { stop("arg result should be a list") }
  if(!(is.R6(spec1) && class(spec1)[1]=="MSSpectrum" && is.R6(spec2) && class(spec2)[1]=="MSSpectrum")) { stop("spec1 and spec2 be MSSpectrum object") }
  result1 <- r_to_py(result)
  private$$py_obj$$getSpectrumAlignment(result1,r_to_py(spec1),r_to_py(spec2))
  tryCatch({
    eval.parent(substitute(result <- py_to_r(result1)))
    invisible()
           }, error = function(e) { invisible()})
})