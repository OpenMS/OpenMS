



    def getData(self):
        """
        Retrieve the data as a NumPy array.

        This method creates a copy of the underlying data, so it's safe to use
        even after the original object is deleted or modified.

        Returns:
            np.ndarray: A NumPy array containing a copy of the data.

        Example:
            >>> da = OSBinaryDataArray()
            >>> da.data = [1.0, 2.0, 3.0]
            >>> arr = da.getData()
            >>> print(arr)  # [1. 2. 3.]
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr.copy())
        return retval

    def getData_mv(self):
        """
        Retrieve the data as a memory view for fast direct access.

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A direct view of the underlying data array.

        Warning:
            The returned memory view refers directly to the underlying data.
            DO NOT store the returned memory view for later use after the
            data array object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> da = OSBinaryDataArray()
            >>> da.data = [1.0, 2.0, 3.0]
            >>> view = da.getData_mv()
            >>> # CORRECT: Use immediately
            >>> total = sum(view)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_view = da.getData_mv()
            >>> # del da  # Memory view now points to invalid memory!
            >>> # print(stored_view[0])  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        return arr
