
InspectInfile$$set("public","getModifications",function(){
  d <- py_call(private$$py_obj$$getModifications)
  k <- lapply( py_to_r(py_builtin$$list(d$$keys())), as.character )
  v <- modify_depth( py_to_r(py_builtin$$list(d$$values())), 2, as.character )
  return(collections::dict(v,k))
})