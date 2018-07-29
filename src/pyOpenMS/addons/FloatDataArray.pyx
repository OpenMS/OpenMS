



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
        cdef libcpp_vector[float].iterator it = fda_.begin()
        cdef unsigned int n = fda_.size()
        cdef np.ndarray[np.float32_t, ndim=1] data
        data = np.zeros( (n,), dtype=np.float32)

        cdef int i = 0
        while it != fda_.end():
            data[i] = deref(it)
            inc(it)
            i += 1

        return data

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

        # cdef float * array_start = <float*>data.data
        cdef int N
        N = len(data)
        # fda_.assign(array_start, array_start + N)
        fda_.reserve(N) # allocate space for incoming data
        for i in range(N):
            fda_.push_back(data[i])


