
DIAScoring$$set("public","dia_by_ion_score",function(spectrum,sequence,charge,bseries_score,yseries_score){
  if(!(is.R6(spectrum) && class(spectrum)[1]=="OSSpectrum")) { stop("arg spectrum wrong type") }
  if(!(is.R6(spectrum) && class(spectrum)[1]=="AASequence")) { stop("arg sequence wrong type") }
  if(!(is_scalar_integer(charge) || (charge %% 1 == 0 )) ) { stop("arg charge wrong type") }
  if(!is_scalar_double(bseries_score)) { stop("arg bseries_score wrong type") }
  if(!is_scalar_double(yseries_score)) { stop("arg yseries_score wrong type") }
  input_bseries_score <- r_to_py(bseries_score)
  input_yseries_score <- r_to_py(yseries_score)
  ans <- private$$py_obj$$dia_by_ion_score(r_to_py(spectrum),r_to_py(sequence),as.integer(charge),input_bseries_score,input_yseries_score)
  tryCatch({
        eval.parent(substitute(bseries_score <- py_to_r(input_bseries_score)))
        eval.parent(substitute(yseries_score <- py_to_r(input_yseries_score)))
        return(ans)
           }, error = function(e) {return(ans)})
} )