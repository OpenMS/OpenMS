
#' DPosition1 Interface
#'
#' @section Constructor:
#' DPosition1$$new()
#' DPosition1$$new(2.2)
#' 
#' @section Extracting Position:
#' \code{dp1 <- DPosition$$new(5.46)}
#' \code{dp1[1]} Get the first value
#' 
#' @name DPosition1
NULL

#' @export
DPosition1 <- R6::R6Class(classname="DPosition1",cloneable = FALSE,
    private = list(py_obj = NA),
    public = list(
    initialize = function(a){
      if(missing(a)){
        private$$py_obj <- Pymod$$DPosition1()
      } else {
        if( "python.builtin.object" %in% class(a) && class_to_wrap(a) == "DPosition1"){
          private$$py_obj <- a
        } else if(is_scalar_double(a)){
          private$$py_obj <- Pymod$$DPosition1(a)
        } else {
          stop("wrong argument a")
        }
      }

    }
    )
)

#' @export
`[.DPosition1` <- function(x,ix){
  stopifnot(R6::is.R6(x))
  if(!ix==1) { stop(paste0("invalid index ",ix)) }
  return(x$$.__enclos_env__$$private$$py_obj[0])
}

#' DPosition2 Interface
#'
#' @section Constructor:
#' DPosition2$$new()
#' DPosition2$$new(2.2,1.1)
#' 
#' @section Extracting Position:
#' \code{dp2 <- DPosition$$new(2.2,1.1)}
#' \code{dp2[1]} Get the first value
#' \code{dp2[2]} Get the second value
#' 
#' @name DPosition2
DPosition2 <- R6::R6Class(classname="DPosition2",cloneable = FALSE,
    private = list(py_obj = NA),

    public = list(
    initialize = function(a,b){
      if(missing(a) && missing(b)){
        private$$py_obj <- Pymod$$DPosition1()
      } else if(missing(b)){
        if("python.builtin.object" %in% class(a) && class_to_wrap(a) == "DPosition2"){ private$$py_obj <- a }
        else { stop("wrong argument provided") }
      } else {
        if (!(is_scalar_double(a) && is_scalar_double(b))) { stop("Both arguments must be double!!") }
        private$$py_obj <- Pymod$$DPosition2(a,b)
      }
    }
    )
)

#' @export
`[.DPosition2` <- function(x,ix){
  stopifnot(R6::is.R6(x))
  if(!ix %in% c(1,2)) { stop(paste0("invalid index ",ix)) }
  return(x$$.__enclos_env__$$private$$py_obj[ix-1])
}
