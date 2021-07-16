
QcMLFile$$set("public","map2csv",function(csv_table,separator){
  if( !(is.environment(csv_table) && identical(parent.env(csv_table), asNamespace("collections")) && strsplit(capture.output(csv_table$$print())," ")[[1]][1] == "dict"
      && all(sapply(csv_table$$keys(),as.character))
      && all(sapply(csv_table$$values(),function(v) is.environment(v) && identical(parent.env(v), asNamespace("collections")) && all(sapply(v$$keys(),as.character)) && all(sapply(v$$values(), as.character))  ))
  )) { stop("arg csv_table wrong type")}
  if(!(is.R6(separator) && class(separator)[1]=="String")) { stop("arg separator wrong type") }
  k <- lapply(csv_table$$keys(),function(k) py_builtin$$bytes(k,'utf-8'))
  v <- lapply(csv_table$$values(),function(v) py_dict(lapply(v$$keys(),function(k) py_builtin$$bytes(k,'utf-8')),lapply(v$$values(),function(v) py_builtin$$bytes(v,'utf-8'))))
  d <- py_dict(k,v)
  res <- private$$py_obj$$map2csv(d,r_to_py(separator))
  return(as.character(res))
})