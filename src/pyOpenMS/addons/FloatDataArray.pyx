



    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, (int, long)), 'arg key wrong type'
    
        cdef long _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <float>_r
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, (int, long)), 'arg key wrong type'
        assert isinstance(value, (int, long, float)), 'arg value wrong type'
        cdef long _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _v = (<float>value)
        deref( self.inst.get() )[(<int>key)] = _v

    def get_data(self):
        """
        Gets the raw data for the float data array

        Example usage: 

          fd = pyopenms.FloatDataArray()
          data = fd.get_data()

        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        cdef unsigned int n = fda_.size()
         
        # Obtain a raw ptr to the beginning of the C++ array
        cdef libcpp_vector[float] * vec_ptr = <libcpp_vector[float]*> fda_
        cdef float * raw_ptr =  address(deref(vec_ptr)[0]) 

        # We use a memory view to get the data from the raw data
        # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef float[:] fda_view = <float[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(fda_view) # numpy array refer to the underlying buffer without copy
        return xarr

    def set_data(self, np.ndarray[float, ndim=1, mode="c"] data not None):
        """
        Sets the raw data for the float data array

        Example usage: 

          fd = pyopenms.FloatDataArray()
          data = numpy.array( [1, 2, 3, 5 ,6] ).astype(numpy.float32)
          fd.set_data(data)

        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        fda_.clear()
        # Note: The np.ndarray[float, ndim=1, mode="c"] assures that we get a C-contiguous numpy array
        # See: https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC

        # We use "assign" to directly to copy the numpy array (stored in
        # data.data) into the std::vector<float> in our FloatDataArray object
        cdef float * array_start = <float*>data.data
        cdef int N
        N = data.size
        fda_.assign(array_start, array_start + N)



