

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

