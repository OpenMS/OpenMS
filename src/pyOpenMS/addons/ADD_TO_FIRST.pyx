from String cimport String as _String
cdef shared_ptr[_String] convString(argument_var)
cdef convOutputString(_String s)



cdef inline convOutputString(_String _r):
    # Generic function to convert OpenMS::String to Python unicode strings
    #
    # This assumes that all data encoded in OpenMS::String is encoded using
    # UTF8 (in practice OpenMS does not make any assumption about encoding and
    # OpenMS::String simply contains a byte sequence, that is implicitely
    # assumed to be ASCII encoded).
    #
    cdef char* c_string = _cast_const_away(<char*> _r.c_str())
    cdef Py_ssize_t length = _r.length()

    try:
        py_result = c_string[:length].decode('UTF-8')
    except UnicodeDecodeError:
        py_result = c_string[:length]

    return py_result

cdef inline shared_ptr[_String] convString(argument_var):
    # Generic function to convert Python strings to OpenMS::String
    # 
    # This allows us to either extract the already existing shared_ptr from a
    # Python-type OpenMS::String holder or create a new shared_ptr from a
    # Python string.
    # In case the user only provides a str, unicode or bytes argument, we will
    # have to create a new OpenMS::String object, the shared_ptr will
    # automatically delete it after it has been used for a function call
    # (presumably the user does not want to have the value returned by
    # reference).
    #
    # Note: even though Python 3.x does not know the unicode keyword,
    # Cython does and uses the PyUnicode_Check call. 
    #
    cdef char* c_string_argument_var
    cdef shared_ptr[_String] res
    if isinstance(argument_var, String):
        res = (<String>argument_var).inst
        return res
    elif isinstance(argument_var, bytes):
        # Simple, convert directly to char* - however there may be zero bytes
        # in the array when using certain encodings, we therefore need to make
        # sure we capture the full length.
        res = shared_ptr[_String](new _String(<char*>argument_var, len(argument_var)))
        return res
    elif isinstance(argument_var, str) or isinstance(argument_var, unicode):
        # First encode with UTF8 (note that output encoding as well as
        # String.toString both decode with UTF8), then convert the result to
        # char*
        py_byte_string = argument_var.encode('UTF-8')
        c_string_argument_var = py_byte_string
        res = shared_ptr[_String](new _String(<char*>c_string_argument_var))
        return res
    else:
        raise Exception("Can only convert the following types to String: pyopenms.String, bytes, str, unicode")

