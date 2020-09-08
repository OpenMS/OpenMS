#' @import reticulate
#' @import R6
#' @import purrr

listDepth <- NULL
Pymod <- NULL
npy <- NULL
py_builtin <- NULL
r_to_py <- NULL


.onLoad <- function(libname, pkgname) {
   reticulate::configure_environment(pkgname)
   Pymod <<- reticulate::import("pyopenms", delay_load = TRUE)
   npy <<- reticulate::import("numpy", convert = F, delay_load = TRUE)
   r_to_py <<- reticulate::r_to_py
   py_builtin <<- reticulate::import_builtins(convert = F)
   listDepth <<- plotrix::listDepth
}

# R6 class object conversion to underlying python object.
`r_to_py.R6` <- function(i,...){
   tryCatch({
       i$.__enclos_env__$private$py_obj
   }, error = function(e) { "conversion not supported for this class"}
   )
}

# Returns the name of wrapper R6 class
class_to_wrap <- function(py_ob){
       class <- tail(strsplit(class(py_ob)[1],"\\.")[[1]],n = 1)
       # To correctly return the class name for Interfaces (BinaryDataArray,Chromatogram,Spectrum) by removing "_Interfaces_"
       comp <- strsplit(class,"_Interfaces_")[[1]]
       if (length(comp) == 1 && comp[1] == class){
           return(class)
       }
       else { return(comp[-1]) }
} 
