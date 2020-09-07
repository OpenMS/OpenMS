

    # This goes into the class
    def _init_0(self):
        self.inst = shared_ptr[_String](new _String())
    
    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==1):
            self.inst = convString(args[0])
        else:
             raise Exception('can not handle type of %s' % (args,)) 

    def toString(self):
        """Cython signature: str toString()
        -- Note: this returns a unicode string and assumes the input is UTF8 encoded.
        """
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

    def c_str(self):
        """Cython signature: const_char * c_str()"""
        # See https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        #    py_string = <bytes> c_string
        # This creates a Python byte string object that holds a copy of the
        # original C string. It can be safely passed around in Python code, and
        # will be garbage collected when the last reference to it goes out of
        # scope. It is important to remember that null bytes in the string act
        # as terminator character, as generally known from C. The above will
        # therefore only work correctly for C strings that do not contain null
        # bytes.
        cdef const_char  * _r = _cast_const_away(self.inst.get().c_str())
        cdef Py_ssize_t length = self.inst.get().length()
        py_result = _r[:length] # This will work correctly also if the char array contains null bytes
        return py_result

