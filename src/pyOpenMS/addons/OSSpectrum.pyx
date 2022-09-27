

    def getMZArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getMZArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def getIntensityArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def getDriftTimeArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL: return None
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getDriftTimeArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL: return None
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def setMZArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setMZArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0)

