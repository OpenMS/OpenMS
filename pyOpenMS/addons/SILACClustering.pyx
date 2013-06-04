

    property rt_min:
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_min
            py_result = <double>_r
            return py_result
    
    property rt_max_spacing:
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_max_spacing
            py_result = <double>_r
            return py_result
