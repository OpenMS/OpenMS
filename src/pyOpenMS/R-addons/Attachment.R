
Attachment$$set("active","tableRows", function(tableRows) {
  if(!missing(tableRows)){
    if( !(is_list(tableRows) && all(sapply(tableRows,function(v0) is.vector(v0) && all(sapply(v0,is_scalar_character)) )) )) { stop("arg tableRows wrong type") }
    private$$py_obj$$tableRows <- map_depth(tableRows,2,function(t) py_builtin$$bytes(t,'utf-8'))
  } else {
        m <- modify_depth(private$$py_obj$$tableRows,2,as.character)
        return(lapply(m, function(v) unlist(v)))
  }
} )

