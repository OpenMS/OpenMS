
Param$$set("public","asDict",function(){
 d <- py_call(private$$py_obj$$asDict)
 k <- lapply(py_to_r(py_builtin$$list(d$$keys())),as.character)
 v <-  py_to_r(py_builtin$$list(d$$values()))
 return(collections::dict(v,k))
})

Param$$set("public","keys",function() {
  res <- lapply(private$$py_obj$$keys(),as.character)
  return(res)
})

Param$$set("public","items",function() {
  res <- purrr::map2(self$$keys(),self$$values(),function(i,j) list(i,j))
  return(res)
})

Param$$set("public","values",function() {
  v <- private$$py_obj$$values()
  return(v)
})

Param$$set("public","update",function(...) {
  arg_list <- list(...)
  if (length(arg_list) == 1){
        if( is.environment(arg_list[[1]]) && identical(parent.env(arg_list[[1]]), asNamespace("collections")) && strsplit(capture.output(arg_list[[1]]$$print())," ")[[1]][1] == "dict") {
          key <- lapply(arg_list[[1]]$$keys(),function(k) py_builtin$$bytes(k,'utf-8'))
          val <- arg_list[[1]]$$values()
          d <- py_dict(key,val)
          private$$py_obj$$update(d)
          invisible()
        } else if( is.R6(arg_list[[1]]) && class(arg_list[[1]]) == "Param" ) {
            private$$py_obj$$update(r_to_py(arg_list[[1]]))
            invisible()
        } else {
          stop("Cannot handle this parameter")
        }
  } else if (length(arg_list)==2) {
    if ( is.R6(arg_list[[1]]) && class(arg_list[[1]])=="Param" || (is.numeric(arg_list[[1]]) && length(arg_list[[1]])==1 && arg_list[[1]] %% 1 == 0)  ) {
      private$$py_obj$$update(r_to_py(arg_list[[1]]),as.integer(arg_list[[2]]))
      invisible()
    } else {
          stop("Cannot handle this parameter")
    }
  } else {
    stop("Invalid parameters provided")
  }
})

Param$$set("public","get",function(key,default) {
  if(!is_scalar_character(key)) { stop("wrong arg key") }
  if(missing(default)) { private$$py_obj$$get(py_builtin$$bytes(key,'utf-8')) }
  else {
    private$$py_obj$$get(py_builtin$$bytes(key,'utf-8'),default)
  }
})

#' @export
`[.Param` <- function(x,ix){
  stopifnot(R6::is.R6(x))
  if(!is_scalar_character(ix)) { stop("key must be a string")}
  res <- x$$.__enclos_env__$$private$$py_obj[py_builtin$$bytes(ix,'utf-8')]
  return(res)
}

#' @export
`[<-.Param` <- function(x,ix,value){
  stopifnot(R6::is.R6(x))
  if(!is_scalar_character(ix)) { stop("key must be a string")}
  x$$.__enclos_env__$$private$$py_obj$$`__setitem__`(py_builtin$$bytes(ix,'utf-8'),value)
  x
}




