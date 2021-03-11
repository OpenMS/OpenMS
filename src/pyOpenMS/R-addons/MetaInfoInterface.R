
MetaInfoInterface$$set("public","getMetaValues",function(){
  d <- py_call(private$$py_obj$$getMetaValues)
  k <- lapply(py_to_r(py_builtin$$list(d$$keys())),as.character)
  v <- py_to_r(py_builtin$$list(d$$values()))
  # v : list(bytes)/list(int)/list(double)/list(list(int))/list(list(double))/list(list(bytes))
  is_nested <- all(sapply(v,function(v1) is_list(v1)))
  if(is_nested){
    v <- map_depth(v,2,function(s){
       if(class(s)[1]=="python.builtin.bytes") as.character(s)
       else { s }
    })
  } else {
    v <- map_depth(v,1, function(s){
       if(class(s)[1]=="python.builtin.bytes") as.character(s)
       else { s }
    })
  }
  return(collections::dict(v,k))
}
)

MetaInfoInterface$$set("public","setMetaValues",function(mmap){
    if( !(is.environment(mmap) && identical(parent.env(mmap), asNamespace("collections")) && strsplit(capture.output(mmap$$print())," ")[[1]][1] == "dict"
      && all(sapply(mmap$$keys(),function(k) is_scalar_character(k)))
      && all(sapply(mmap$$values(),function(v) (is_scalar_integer(v) || is_scalar_double(v) || is_scalar_character(v)) || (is_list(v) && all(lapply(v,function(vi) (is_scalar_integer(vi) || is_scalar_double(vi) || is_scalar_character(vi))))) )) )
  ) { stop("arg mmap wrong type")}
  # handle corner case of list(list(char)) conversion to list(list(bytes)) as python expects a list of bytes for StringList
  is_nested_char <- all(sapply(mmap$$values(),function(v) is_list(v) && all(sapply(v,function(vi) is_scalar_character(vi))) ))
  v <- mmap$$values()
  k <- mmap$$keys()
  if(is_nested_char) {
    v <- map_depth(v,2,function(vi) py_builtin$$bytes(vi,'utf-8'))
  }
  private$$py_obj$$setMetaValues(py_dict(k,v))
  invisible()
})