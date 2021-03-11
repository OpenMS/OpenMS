
ConsensusMap$$set("public","setUniqueIds",function() {
  private$$py_obj$$setUniqueIds()
} )

ConsensusMap$$set("public","getColumnHeaders",function() {
  d <- py_call(private$$py_obj$$getColumnHeaders)
  key <- py_to_r(py_builtin$$list(d$$keys()))
  val <- py_to_r(py_builtin$$list(d$$values()))
  return (collections::dict(lapply(val,function(v) ColumnHeader$$new(v)),key))
} )

ConsensusMap$$set("public","setColumnHeaders",function(in_0) {
    if( !(is.environment(in_0) && identical(parent.env(in_0), asNamespace("collections")) && strsplit(capture.output(in_0$$print())," ")[[1]][1] == "dict"
      && all(sapply(map$$keys(),function(k) is_scalar_integer(k) || k %% 1 == 0 ))
      && all(sapply(map$$values(),function(v) is.R6(v) && class(v)[1]=="ColumnHeader" )))
  ) { stop("arg map wrong type")}
  d <- py_dict(lapply(in_0$$keys(),as.integer),lapply(in_0$$values(),function(v) {r_to_py(v)}))
  private$$py_obj$$setColumnHeaders(d)
  k <- py_to_r(py_builtin$$list(d$$keys()))
  v <- lapply(py_to_r(py_builtin$$list(d$$values())),function(c) ColumnHeader$$new(c))
  tryCatch({
             eval.parent(substitute(in_0 <- collections::dict(v,k)))
             invisible()
        }, error = function(c) { invisible() }
  )

} )
