
MSQuantifications$$set("public","registerExperiment", function(exp,labels) {
  if(!(is.R6(exp) && class(exp)[1]=="MSExperiment")) { stop("arg exp wrong type") }
  if (!(is_list(labels) && all(sapply(labels), function(l2) is_list(l2) && all(sapply(l2, function(l3) is_list(l3) && all(sapply(l3, function(l4) is_list(l4) && length(l4)==2 && all(class(l4[[1]])==c("String","R6")) && is_scalar_double(l4[[2]]))))))) ) { stop("arg labels wrong type") }
  private$$py_obj$$registerExperiment(exp,modify_depth(labels, 4,function(a) r_to_py(a)))
  invisible()
} )