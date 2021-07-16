
#' @export
`[.StringDataArray` <- function(x,ix){
  stopifnot(is.R6(x))
  if(!isTRUE(all.equal(ix,as.integer(ix)))) { stop("index should be an integer") }
  tryCatch({
    res <- x$$.__enclos_env__$$private$$py_obj[as.integer(ix-1)]
    return(as.character(res))
           }, error = function(e){ e })
}

#' @export
`[<-.StringDataArray` <- function(x,ix,value) {
  stopifnot(is.R6(x))
  if(!isTRUE(all.equal(ix,as.integer(ix)))) { stop("index should be an integer") }
  if(!( (is.R6(value) && class(value)[1]=="String") || is_scalar_character(value) )) { stop("value must be a string or String object") }
  tryCatch({
    x$$.__enclos_env__$$private$$py_obj$$`__setitem__`(as.integer(ix-1),r_to_py(value))
    invisible(x)
           }, error = function(e){ e })
}