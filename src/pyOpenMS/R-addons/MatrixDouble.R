
# Returns a matrix corresponding to 2d numpy array
MatrixDouble$$set("public","get_matrix",function(){
  private$$py_obj$$get_matrix()
} )

# Returns a matrix corresponding to 2d numpy array
MatrixDouble$$set("public","get_matrix_as_view",function(){
  private$$py_obj$$get_matrix_as_view()
} )

MatrixDouble$$set("public","set_matrix",function(data){
  if(!(is.matrix(data) && all(sapply(seq_along(NROW(data)), function(d) is_double(data[d,]) ))) ) { stop("arg data wrong type") }
  private$$py_obj$$set_matrix(npy$$array(data, order = "c"))
  invisible()
} )
