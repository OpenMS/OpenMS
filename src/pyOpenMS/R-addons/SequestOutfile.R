
SequestOutfile$$set("public","getSequences",function(database_filename,ac_position_map,sequences,found,not_found){
  if(!(is.R6(database_filename) && class(database_filename)[1]=="String")) { stop("arg database_filename wrong type") }
  if( !(is.environment(ac_position_map) && identical(parent.env(ac_position_map), asNamespace("collections")) && strsplit(capture.output(ac_position_map$$print())," ")[[1]][1] == "dict"
      && all(sapply(ac_position_map$$keys(),is_scalar_character))
      && all(sapply(ac_position_map$$values(),function(v) isTRUE(all.equal(v,as.integer(v))) ))
  )) { stop("arg ac_position_map wrong type")}
  if(!( (is_list(sequences) || is_vector(sequences)) && all(sapply(sequences,is_scalar_character)) )) { stop("arg sequences wrong type") }
  if(!( is_list(found) && all(sapply(found,function(f) is_list(f,2) && is_scalar_character(f[[1]]) && is_scalar_integer(f[[2]]))) )) { stop("arg found wrong type") }
  if( !(is.environment(not_found) && identical(parent.env(not_found), asNamespace("collections")) && strsplit(capture.output(not_found$$print())," ")[[1]][1] == "dict"
      && all(sapply(not_found$$keys(),is_scalar_character))
      && all(sapply(not_found$$values(),function(v) isTRUE(all.equal(v,as.integer(v))) ))
  )) { stop("arg not_found wrong type") }
  
  k1 <- lapply(ac_position_map$$keys(),function(k) py_builtin$$bytes(k,'utf-8'))
  v1 <- lapply(ac_position_map$$values(),as.integer)
  temp1 <- py_dict(k1,v1)
  
  seq1 <- r_to_py(lapply(sequences,function(s) py_builtin$$bytes(s,'utf-8')))
  found1 <- r_to_py(modify_depth(found,1,function(p) list(py_builtin$$bytes(p[[1]],'utf-8'),p[[2]])))
  
  k1 <- lapply(not_found$$keys(),function(k) py_builtin$$bytes(k,'utf-8'))
  v1 <- lapply(not_found$$values(),as.integer)
  temp2 <- py_dict(k1,v1)
  
  private$$py_obj$$getSequences(database_filename,temp1,seq1,found1,temp2)
  k1 <- lapply(py_to_r(py_builtin$$list(temp1$$keys())),as.character)
  v1 <- py_to_r(py_builtin$$list(temp1$$values()))
  temp1 <- collections::dict(v1,k1)
  
  seq1 <- lapply(py_to_r(seq1),as.character)
  found1 <- modify_depth(py_to_r(found1),1,function(p) list(as.character(p[[1]]),p[[2]]))
  k1 <- lapply(py_to_r(py_builtin$$list(temp2$$keys())),as.character)
  v1 <- py_to_r(py_builtin$$list(temp2$$values()))
  temp2 <- collections::dict(v1,k1)
  
  tryCatch({
    eval.parent(ac_position_map <- temp1)
    eval.parent(sequences <- seq1)
    eval.parent(found <- found1)
    eval.parent(not_found <- temp2)
    invisible()
           }, error = function(){ invisible() })

})