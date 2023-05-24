

    def getMZArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef np.ndarray[np.float64_t, ndim=1] retval
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr.copy())
        return retval

    def getMZArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def getIntensityArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef np.ndarray[np.float64_t, ndim=1] retval
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr.copy())
        return retval

    def getIntensityArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def getDriftTimeArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL: return None
        cdef np.ndarray[np.float64_t, ndim=1] retval
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr.copy())
        return retval

    def getDriftTimeArray_mv(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL: return None
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr

    def getDataArrays(self):
        cdef list py_result = []
        cdef OSBinaryDataArray rv

        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ]  _r = self.inst.get().getDataArrays()

        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ].iterator it = _r.begin()
        while it != _r.end():
            rv = OSBinaryDataArray.__new__(OSBinaryDataArray)
            rv.inst = deref(it)
            py_result.append(rv)
            inc(it)
        return py_result

    def setDataArrays(self, list inp):
        assert isinstance(inp, list) and all(isinstance(ele_inp, OSBinaryDataArray) for ele_inp in inp), 'Input has to be a list of elements of type OSBinaryDataArray'


        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ]  _r = self.inst.get().getDataArrays() 
        _r.clear() 

        cdef OSBinaryDataArray rv
        for rv in inp:
            _r.push_back(rv.inst)

        self.inst.get().setDataArrays(_r)

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

