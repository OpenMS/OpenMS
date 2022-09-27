

    def getTimeArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getDataArrays(self):
        cdef list py_result = []
        cdef OSBinaryDataArray rv

        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ]  _r = self.inst.get().getDataArrays()
        cdef shared_ptr[_OSBinaryDataArray] v0

        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ].iterator it = _r.begin()
        while it != _r.end():
            v0 = deref(it)
            rv = OSBinaryDataArray.__new__(OSBinaryDataArray)
            rv.inst = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray(deref(v0.get())))
            py_result.append(rv)
            inc(it)
        return py_result

    def setTimeArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setTimeArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0)

