

    def getTimeArray(self):
        """
        getTimeArray(self: Interfaces_Chromatogram) -> List[float]
        
        Get the time array as a list.
        """
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        """
        getIntensityArray(self: Interfaces_Chromatogram) -> List[float]
        
        Get the intensity array as a list.
        """
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def setTimeArray(self, list data):
        """
        setTimeArray(self: Interfaces_Chromatogram, data: List[float]) -> None
        
        Set the time array from a list.
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setTimeArray(v0)

    def setIntensityArray(self, list data):
        """
        setIntensityArray(self: Interfaces_Chromatogram, data: List[float]) -> None
        
        Set the intensity array from a list.
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0)

