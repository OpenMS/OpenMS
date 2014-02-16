


    cdef shared_ptr[_IndexedMzMLDecoder] inst
 
    def __dealloc__(self):
        self.inst.reset()

    def findIndexListOffset(self, bytes in_ ,  buffersize ):
        assert isinstance(in_, bytes), 'arg in_ wrong type'
        assert isinstance(buffersize, (int, long)), 'arg buffersize wrong type'
    
        cdef long long _r = self.inst.get().findIndexListOffset((_String(<char *>in_)), (<int>buffersize))
        py_result = <long long>_r
        return py_result 
