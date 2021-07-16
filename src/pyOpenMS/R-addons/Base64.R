
# Cython signature: void encode64(libcpp_vector[double] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)
Base64$$set("public","encode64",function(in_,to_byte_order,out,zlib_compression) {
  if (!(is.vector(in_) && all(sapply(in_, function(i) is_scalar_double(i))))) { stop("arg in_ wrong type") }
  if (!to_byte_order %in% c(0,1)) { stop("arg to_byte_order wrong type") }
  if (!(is.R6(out) && class(out)[1] == "String")) { stop("arg out wrong type") }
  if (!is.logical(zlib_compression)) { stop("arg zlib_compression wrong type") }
  private$$py_obj$$encode64(as.list(in_),as.integer(to_byte_order),r_to_py(out),as.integer(zlib_compression))
} )

# Cython signature: void decode64(const String & in_, ByteOrder from_byte_order, libcpp_vector[double] & out, bool zlib_compression)
Base64$$set("public","decode64",function(in_,from_byte_order,out,zlib_compression) {
  if ( !(is_scalar_character(in_) || is.R6(in_) && class(in_)[1] == "String") ) { stop("arg in_ wrong type") }
  if (!from_byte_order %in% c(0,1)) { stop("arg from_byte_order wrong type") }
  if (!(is.vector(out) && all(sapply(out, function(o) is_scalar_double(o))) )) { stop("arg out wrong type") }
  if (!is.logical(zlib_compression)) { stop("arg zlib_compression wrong type") }
  out1 <- r_to_py(as.list(out))
  private$$py_obj$$decode64(r_to_py(in_),as.integer(from_byte_order),out1,as.integer(zlib_compression))

  tryCatch({
             eval.parent(substitute(out <- py_to_r(out1)))
             invisible()
        }, error = function(c) { invisible() }
  )
} )

# Cython signature: void encode32(libcpp_vector[float] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)
Base64$$set("public","encode32",function(in_,to_byte_order,out,zlib_compression) {
  if (!(is.vector(in_) && all(sapply(in_, function(i) is_scalar_double(i))))) { stop("arg in_ wrong type") }
  if (!to_byte_order %in% c(0,1)) { stop("arg to_byte_order wrong type") }
  if (!(is.R6(out) && class(out)[1] == "String")) { stop("arg out wrong type") }
  if (!is.logical(zlib_compression)) { stop("arg zlib_compression wrong type") }
  private$$py_obj$$encode32(as.list(in_),as.integer(to_byte_order),r_to_py(out),as.integer(zlib_compression))
} )

# Cython signature: void decode32(const String & in_, ByteOrder from_byte_order, libcpp_vector[float] & out, bool zlib_compression)
Base64$$set("public","decode32",function(in_,from_byte_order,out,zlib_compression) {
  if ( !(is_scalar_character(in_) || is.R6(in_) && class(in_)[1] == "String") ) { stop("arg in_ wrong type") }
  if (!from_byte_order %in% c(0,1)) { stop("arg from_byte_order wrong type") }
  if (!(is.vector(out) && all(sapply(out, function(o) is_scalar_double(o))) )) { stop("arg out wrong type") }
  if (!is.logical(zlib_compression)) { stop("arg zlib_compression wrong type") }
  out1 <- r_to_py(as.list(out))
  private$$py_obj$$decode32(r_to_py(in_),as.integer(from_byte_order),out1,as.integer(zlib_compression))

  tryCatch({
             eval.parent(substitute(out <- py_to_r(out1)))
             invisible()
        }, error = function(c) { invisible() }
  )
} )










