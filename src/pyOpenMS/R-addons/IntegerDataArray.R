
# Gets the raw data for the float data array
# Example usage:
# fd = IntegerDataArray$$new()
# data = fd$$get_data()
IntegerDataArray$$set("public","get_data",function(){
  ans <- private$$py_obj$$get_data()
  return(as.vector(ans))
}
)


# Sets the raw data for the float data array
#
# Example usage:
#
# fd = IntegerDataArray$$new()
# data = 1:5
# fd$$set_data(data)
IntegerDataArray$$set("public","set_data",function(data){
  if (!( is_vector(data) && all(sapply(data, function(d) isTRUE(all.equal(d,as.integer(d))))) && (is.null(ncol(data)) || is.na(ncol(data))) )) { stop(paste0("Wrong argument ",data)) }
  private$$py_obj$$set_data(npy$$array(data)$$astype(npy$$intc))
}
)

#' @export
`[.IntegerDataArray` <- function(x,ix){
  stopifnot(R6::is.R6(x))
  if(!(ix %% 1 == 0)) { stop("index must be integer") }
  tryCatch({
    return(x$$.__enclos_env__$$private$$py_obj[as.integer(ix-1)])
           }, error = function(e) { stop(paste0("invalid index ",ix)) }
  )
}

#' @export
`[<-.IntegerDataArray` <- function(x,ix,value){
  stopifnot(R6::is.R6(x))
  if(!(ix %% 1 == 0)) { stop("index must be integer") }
  if(!(value %% 1 == 0)) { stop("value must be integer") }
    tryCatch({
    x$$.__enclos_env__$$private$$py_obj$$`__setitem__`(as.integer(ix-1), as.integer(value))
    invisible(x)
           }, error = function(e) { stop(paste0("invalid index ",ix)) }
  )
}

