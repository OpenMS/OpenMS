



    def __getitem__(self,  in_0 ):
        """
        __getitem__(self: FloatDataArray, in_0: int) -> float
        
        Get element at index in_0.
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert in_0 >= 0, 'arg in_0 cannot be negative'

        cdef unsigned int _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <float>_r
        return py_result

    def __setitem__(self, key, value):
        """
        __setitem__(self: FloatDataArray, key: int, value: Union[int, float]) -> None
        
        Set element at index key to value.
        """
        assert isinstance(key, int), 'arg key wrong type'
        assert isinstance(value, (int, float)), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _v = (<float>value)
        deref( self.inst.get() )[(<int>key)] = _v

    def get_data(self):
        """
        get_data(self: FloatDataArray) -> np.ndarray
        
        Gets a copy of the data as a numpy array (safe).

        This method creates a copy of the underlying data, so it's safe to use
        even after the original FloatDataArray object is deleted or modified.

        Returns:
            np.ndarray: A numpy array (float32) containing a copy of the data.

        Example usage:

        .. code-block:: python

            fd = pyopenms.FloatDataArray()
            fd.push_back(1.0)
            fd.push_back(2.0)
            data = fd.get_data()  # Safe copy
            del fd  # Original can be deleted, data is still valid
            print(data)  # [1.0, 2.0]
        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        cdef unsigned int n = fda_.size()

        if n == 0:
            return np.empty(0, dtype=np.float32)

        cdef np.ndarray[np.float32_t, ndim=1] result = np.empty(n, dtype=np.float32)
        cdef unsigned int i
        for i in range(n):
            result[i] = deref(fda_)[i]

        return result

    def get_data_mv(self):
        """
        get_data_mv(self: FloatDataArray) -> Optional[np.ndarray]
        
        Gets the raw data for the float data array as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        .. warning::
           This returns a fast but unsafe view on the underlying vector data.
           Make sure that the object you are getting the data from survives the
           usage of this view! Specifically, do NOT use it like this:
           :code:`data = spectrum.getFloatDataArrays()[0].get_data_mv()` since the
           underlying FloatDataArray is temporary for this line and will most
           likely be garbage collected right after that line.

           For safe access, use get_data() instead.

        Returns:
            np.ndarray: A numpy array that is a view of the underlying data.
                       Returns None if no data present.

        Example usage:

        .. code-block:: python

            fd = pyopenms.FloatDataArray()
            fd.push_back(1.0)
            fd.push_back(2.0)
            data = fd.get_data_mv()  # Memory view - changes affect original
        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        cdef unsigned int n = fda_.size()

        if n == 0:
            return None

        # Obtain a raw ptr to the beginning of the C++ array
        cdef libcpp_vector[float] * vec_ptr = <libcpp_vector[float]*> fda_
        cdef float * raw_ptr = address(deref(vec_ptr)[0])

        # We use a memory view to get the data from the raw data
        # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html
        # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef float[:] fda_view = <float[:n]>raw_ptr  # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(fda_view)  # numpy array refer to the underlying buffer without copy
        return xarr

        ## # Alternatively, we could also use the assignment method:
        ## cdef _FloatDataArray * fda_ = self.inst.get()
        ## cdef libcpp_vector[float].iterator it = fda_.begin()
        ## cdef unsigned int n = fda_.size()
        ## cdef np.ndarray[np.float32_t, ndim=1] data
        ## data = np.zeros( (n,), dtype=np.float32)
        ## 
        ## cdef int i = 0
        ## while it != fda_.end():
        ##     data[i] = deref(it)
        ##     inc(it)
        ##     i += 1
        ## 
        ## return data
        ##

        # We see an improvement of ca 50x using the memory view method (and a
        # 2600x fold improvement compared to pure Python).

    def set_data(self, np.ndarray[float, ndim=1, mode="c"] data not None):
        """
        set_data(self: FloatDataArray, data: np.ndarray) -> None
        
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


        ## # Alternatively, we could also "push_back" each element:
        ## fda_.reserve(N)
        ## for i in range(N):
        ##     fda_.push_back(data[i])
        
        # We see an improvement of ca 10x using the assign method (and an 600x
        # improvement compared to pure Python).

    def __len__(self):
        """
        __len__(self: FloatDataArray) -> int
        
        Return the number of elements in the array.
        """
        return self.inst.get().size()

    def __str__(self):
        """
        __str__(self: FloatDataArray) -> str
        
        Return a string representation of the FloatDataArray object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        __repr__(self: FloatDataArray) -> str
        
        Return a string representation of the FloatDataArray object.

        Returns key properties in a readable format:
        FloatDataArray(name='ion_mobility', size=100)
        """
        cdef unsigned int n = self.inst.get().size()
        
        parts = []
        
        # Get name from MetaInfoDescription
        name = self.getName()
        if name:
            parts.append(f"name='{name}'")
        
        parts.append(f"size={n}")
        
        return f"FloatDataArray({', '.join(parts)})"
