from String cimport String as _String
cdef class StringPtrWrapper:
    cdef shared_ptr[_String] wrapped
    def __dealloc__(self): self.wrapped.reset()
cdef convString(argument_var):
    cdef StringPtrWrapper wrapper = StringPtrWrapper()
    #cdef _String * cs_argument_var
    cdef char* c_string_argument_var
    # cdef shared_ptr[_String] res
    if isinstance(argument_var, String):
        wrapper.wrapped = (<String>argument_var).inst
        return wrapper
    elif isinstance(argument_var, bytes):
        wrapper.wrapped = shared_ptr[_String](new _String(<char*>argument_var))
        return wrapper
    elif isinstance(argument_var, str):
        py_byte_string = argument_var.encode('UTF-8')
        c_string_argument_var = py_byte_string
        # cs_argument_var = new _String((<char *>c_string_argument_var))
        wrapper.wrapped = shared_ptr[_String](new _String(<char*>c_string_argument_var))
        return wrapper
    elif isinstance(argument_var, unicode):
        py_byte_string = argument_var.encode('UTF-8')
        c_string_argument_var = py_byte_string
        #cs_argument_var = new _String((<char *>c_string_argument_var))
        wrapper.wrapped = shared_ptr[_String](new _String(<char*>c_string_argument_var))
        return wrapper
    else:
        raise Exception("Can only convert the following types to String: pyopenms.String, bytes, str, unicode")



    def _init_0(self):
        self.inst = shared_ptr[_String](new _String())
    
    def _init_1(self, String in_0 ):
        assert isinstance(in_0, String), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_String](new _String(deref(in_0.inst.get())))

    def _init_2(self, bytes in_0 ):
        assert isinstance(in_0, bytes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_String](new _String((<char *>in_0)))

    def _init_3(self, str in_0 ):
        assert isinstance(in_0, str), 'arg in_0 wrong type'
    
        py_byte_string = in_0.encode('UTF-8')
        cdef char* c_string = py_byte_string
        self.inst = shared_ptr[_String](new _String((<char *>c_string)))

    def _init_4(self, unicode in_0 ):
        assert isinstance(in_0, unicode), 'arg in_0 wrong type'
    
        py_byte_string = in_0.encode('UTF-8')
        cdef char* c_string = py_byte_string
        self.inst = shared_ptr[_String](new _String((<char *>c_string)))
    
    def __init__(self, *args):

        #cdef shared_ptr[_String] mystr = (<shared_ptr[_String]>convString(*args))
        cdef StringPtrWrapper mystr = (<StringPtrWrapper>convString(args[0]))
        cdef shared_ptr[_String] mystr2 = (<StringPtrWrapper>convString(args[0])).wrapped


        # Note: even though Python 3.x does not know the unicode keyword,
        # Cython does and uses the PyUnicode_Check call. We just need to check
        # for unicode last, this will ensure all Python 3.x call go through the
        # "str" branch below:
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], String)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], bytes)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], str)):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], unicode)):
             self._init_4(*args)
        else:
             raise Exception('can not handle type of %s' % (args,)) 

    def toString(self):
        # Decodes the C string to unicode
        cdef char* c_string = _cast_const_away(self.inst.get().c_str())
        cdef Py_ssize_t length = self.inst.get().length()
        ustring = c_string[:length].decode('UTF-8')
        return ustring

    def __bytes__(self):
        return self.c_str()

    def __str__(self):
        # In Python 2 this is supposed to return bytes, in Python 3 it is
        # supposed to return Unicode. We cannot really support both, so we
        # return bytes (as a str object for Python 3).
        # Dont use this in serious code.
        return str(self.c_str())

    def __unicode__(self):
        return self.toString()

    def __repr__(self):
        return self.c_str()

    def __hash__(self):
        return hash(self.c_str())

