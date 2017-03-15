



    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, (int, long)), 'arg in_0 wrong type'
    
        cdef long _idx = (<int>in_0)
        cdef int _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <int>_r
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, (int, long)), 'arg key wrong type'
        assert isinstance(value, (int, long)), 'arg value wrong type'
        cdef long _idx = (<int>key)
        cdef int _v = (<int>value)
        deref(self.inst.get())[(<int>key)] = _v

