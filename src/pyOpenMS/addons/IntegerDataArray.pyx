



    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, (int, long)), 'arg in_0 wrong type'
    
        cdef long _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <int>_r
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, (int, long)), 'arg key wrong type'
        assert isinstance(value, (int, long)), 'arg value wrong type'
        cdef long _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _v = (<int>value)
        deref(self.inst.get())[(<int>key)] = _v

    def get_data(self):
        """
        Gets the raw data for the integer data array

        Example usage: 

          idata = pyopenms.IntegerDataArray()
          data = idata.get_data()

        """

        cdef _IntegerDataArray * ida_ = self.inst.get()
        cdef libcpp_vector[int].iterator it = ida_.begin()
        cdef unsigned int n = ida_.size()
        cdef np.ndarray[int, ndim=1] data
        data = np.zeros( (n,), dtype=np.intc)

        cdef int i = 0
        while it != ida_.end():
            data[i] = deref(it)
            inc(it)
            i += 1

        return data

    def set_data(self, np.ndarray[int, ndim=1, mode="c"] data not None):
        """
        Sets the raw data for the integer data array

        Example usage: 

          idata = pyopenms.IntegerDataArray()
          data = numpy.array( [1, 2, 3, 5 ,6] ).astype(np.intc)
          idata.set_data(data)

        """

        cdef _IntegerDataArray * ida_ = self.inst.get()
        ida_.clear()

        cdef int N
        N = len(data)
        ida_.reserve(N) # allocate space for incoming data
        for i in range(N):
            ida_.push_back(data[i])

