

#' @title Streampos object to integer.
#' @export
`as.integer.streampos` <- function(x, ...) {
    py_to_r(py_builtin$$int(r_to_py(x)))
}