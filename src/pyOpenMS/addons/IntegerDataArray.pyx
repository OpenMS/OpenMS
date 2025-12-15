



    def __getitem__(self,  in_0 ):
        """
        __getitem__(self: IntegerDataArray, in_0: int) -> int
        
        Get element at index in_0.
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert in_0 >= 0, 'arg in_0 cannot be negative'

        cdef unsigned int _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <int>_r
        return py_result

    def __setitem__(self, key, value):
        """
        __setitem__(self: IntegerDataArray, key: int, value: int) -> None
        
        Set element at index key to value.
        """
        assert isinstance(key, int), 'arg key wrong type'
        assert isinstance(value, int), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _v = (<int>value)
        deref(self.inst.get())[(<int>key)] = _v

    def get_data(self):
        """
        get_data(self: IntegerDataArray) -> np.ndarray
        
        Gets a copy of the data as a numpy array (safe).

        This method creates a copy of the underlying data, so it's safe to use
        even after the original IntegerDataArray object is deleted or modified.

        Returns:
            np.ndarray: A numpy array (int32) containing a copy of the data.

        Example usage:

        .. code-block:: python

            ida = pyopenms.IntegerDataArray()
            ida.push_back(1)
            ida.push_back(2)
            data = ida.get_data()  # Safe copy
            del ida  # Original can be deleted, data is still valid
            print(data)  # [1, 2]
        """
        cdef _IntegerDataArray * ida_ = self.inst.get()
        cdef unsigned int n = ida_.size()

        if n == 0:
            return np.empty(0, dtype=np.int32)

        cdef np.ndarray[np.int32_t, ndim=1] result = np.empty(n, dtype=np.int32)
        cdef unsigned int i
        for i in range(n):
            result[i] = deref(ida_)[i]

        return result

    def get_data_mv(self):
        """
        get_data_mv(self: IntegerDataArray) -> Optional[np.ndarray]
        
        Gets the raw data for the integer data array as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        .. warning::
           This returns a fast but unsafe view on the underlying vector data.
           Make sure that the object you are getting the data from survives the
           usage of this view! Specifically, do NOT use it like this:
           :code:`data = spectrum.getIntegerDataArrays()[0].get_data_mv()` since the
           underlying IntegerDataArray is temporary for this line and will most
           likely be garbage collected right after that line.

           For safe access, use get_data() instead.

        Returns:
            np.ndarray: A numpy array that is a view of the underlying data.
                       Returns None if no data present.

        Example usage:

        .. code-block:: python

            ida = pyopenms.IntegerDataArray()
            ida.push_back(1)
            ida.push_back(2)
            data = ida.get_data_mv()  # Memory view - changes affect original
        """
        cdef _IntegerDataArray * ida_ = self.inst.get()
        cdef unsigned int n = ida_.size()

        if n == 0:
            return None

        # Obtain a raw ptr to the beginning of the C++ array
        cdef libcpp_vector[int] * vec_ptr = <libcpp_vector[int]*> ida_
        cdef int * raw_ptr = address(deref(vec_ptr)[0])

        # We use a memory view to get the data from the raw data
        # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html
        cdef int[:] ida_view = <int[:n]>raw_ptr  # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(ida_view)  # numpy array refer to the underlying buffer without copy
        return xarr

    def set_data(self, np.ndarray[int, ndim=1, mode="c"] data not None):
        """
        set_data(self: IntegerDataArray, data: np.ndarray) -> None
        
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

    def __len__(self):
        """
        __len__(self: IntegerDataArray) -> int
        
        Return the number of elements in the array.
        """
        return self.inst.get().size()

