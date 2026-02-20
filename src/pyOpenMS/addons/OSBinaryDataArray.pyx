



    def get_data(self):
        """
        Retrieve the data as a NumPy array (copy).

        This method creates a copy of the underlying data, so it's safe to use
        even after the original object is deleted or modified.

        Returns:
            np.ndarray: A NumPy array (float64) containing a copy of the data.

        Example:
            >>> da = OSBinaryDataArray()
            >>> da.data = [1.0, 2.0, 3.0]
            >>> arr = da.get_data()
            >>> print(arr)  # [1. 2. 3.]
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_data_mv(self):
        """
        Retrieve the data as a memory view for fast direct access (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A direct view of the underlying data array, or None if empty.

        Warning:
            The returned memory view refers directly to the underlying data.
            DO NOT store the returned memory view for later use after the
            data array object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> da = OSBinaryDataArray()
            >>> da.data = [1.0, 2.0, 3.0]
            >>> view = da.get_data_mv()
            >>> # CORRECT: Use immediately
            >>> total = sum(view)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_view = da.get_data_mv()
            >>> # del da  # Memory view now points to invalid memory!
            >>> # print(stored_view[0])  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr
