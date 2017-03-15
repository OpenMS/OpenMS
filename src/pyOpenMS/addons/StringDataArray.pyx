



    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, (int, long)), 'arg in_0 wrong type'
    
        cdef long _idx = (<int>in_0)
        cdef _String _r = deref(self.inst.get())[(<int>in_0)]
        py_result = _cast_const_away(<char*>_r.c_str())
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, (int, long)), 'arg key wrong type'
        assert (isinstance(value, str) or isinstance(value, unicode) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'

        cdef long _idx = (<int>key)
        cdef shared_ptr[_String] _s = convString(value)
        deref(self.inst.get())[(<int>key)] = deref(_s.get())

