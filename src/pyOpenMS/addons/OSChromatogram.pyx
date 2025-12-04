

    def get_time_array(self):
        """
        Get the time array of the chromatogram as a numpy array (copy).

        Returns:
            np.ndarray: A 1D numpy array (float64) containing the time values (retention times)
                       for each data point in the chromatogram.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> times = chromatogram.get_time_array()
            >>> print(f"Retention times: {times}")
            >>> # Safe to modify the returned array
            >>> times[0] = 0.0  # This won't affect the original data
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_time_array_mv(self):
        """
        Get the time array of the chromatogram as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A memory view of the time values (retention times)
                       for each data point in the chromatogram, or None if empty.

        Warning:
            DO NOT store the returned memory view for later use after the
            chromatogram object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> times_mv = chromatogram.get_time_array_mv()
            >>> # CORRECT: Use immediately
            >>> max_time = max(times_mv)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_mv = chromatogram.get_time_array_mv()
            >>> # del chromatogram  # Memory view now points to invalid memory!
            >>> # print(stored_mv[0])  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr

    def get_intensity_array(self):
        """
        Get the intensity array of the chromatogram as a numpy array (copy).

        Returns:
            np.ndarray: A 1D numpy array (float64) containing the intensity values
                       for each data point in the chromatogram.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> intensities = chromatogram.get_intensity_array()
            >>> print(f"Max intensity: {max(intensities)}")
            >>> # Safe to modify the returned array
            >>> intensities *= 2.0  # This won't affect the original data
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_intensity_array_mv(self):
        """
        Get the intensity array of the chromatogram as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A memory view of the intensity values
                       for each data point in the chromatogram, or None if empty.

        Warning:
            DO NOT store the returned memory view for later use after the
            chromatogram object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> intensities_mv = chromatogram.get_intensity_array_mv()
            >>> # CORRECT: Use immediately for calculations
            >>> total_intensity = sum(intensities_mv)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_mv = chromatogram.get_intensity_array_mv()
            >>> # del chromatogram  # Memory view now points to invalid memory!
            >>> # print(sum(stored_mv))  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr

    def get_data_arrays(self):
        """
        Get all data arrays associated with the chromatogram.

        This includes the time array, intensity array, and any additional
        meta data arrays that may be present in the chromatogram.

        Returns:
            list: A list of OSBinaryDataArray objects containing all data arrays
                 associated with this chromatogram.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> data_arrays = chromatogram.get_data_arrays()
            >>> print(f"Number of data arrays: {len(data_arrays)}")
            >>> for i, array in enumerate(data_arrays):
            ...     print(f"Array {i}: {len(array.get_data())} data points")
        """
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

    def set_data_arrays(self, list inp):
        """
        Set all data arrays for the chromatogram.

        This method replaces all existing data arrays with the provided ones.
        The input should include time array, intensity array, and any additional
        meta data arrays.

        Args:
            inp (list): A list of OSBinaryDataArray objects to set as the
                       chromatogram's data arrays.

        Raises:
            AssertionError: If inp is not a list or contains elements that are
                           not OSBinaryDataArray instances.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> time_array = OSBinaryDataArray()
            >>> intensity_array = OSBinaryDataArray()
            >>> # ... populate arrays with data ...
            >>> chromatogram.set_data_arrays([time_array, intensity_array])
        """
        assert isinstance(inp, list) and all(isinstance(ele_inp, OSBinaryDataArray) for ele_inp in inp), 'Input has to be a list of elements of type OSBinaryDataArray'

        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ]  _r = self.inst.get().getDataArrays()
        _r.clear()

        cdef OSBinaryDataArray rv
        for rv in inp:
            _r.push_back(rv.inst)

        self.inst.get().setDataArrays(_r)

    def set_time_array(self, list data):
        """
        Set the time array (retention times) for the chromatogram.

        Args:
            data (list): A list of numeric values representing the retention times
                        for each data point in the chromatogram. Values should be
                        in ascending order.

        Raises:
            AssertionError: If data is not a list.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> retention_times = [0.5, 1.0, 1.5, 2.0, 2.5]  # in minutes
            >>> chromatogram.set_time_array(retention_times)
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() )
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setTimeArray(v0)

    def set_intensity_array(self, list data):
        """
        Set the intensity array for the chromatogram.

        Args:
            data (list): A list of numeric values representing the intensity values
                        for each data point in the chromatogram. The length should
                        match the time array length.

        Raises:
            AssertionError: If data is not a list.

        Example:
            >>> chromatogram = OSChromatogram()
            >>> intensities = [1000.0, 1500.0, 2000.0, 1200.0, 800.0]
            >>> chromatogram.set_intensity_array(intensities)
            >>> # Make sure length matches time array
            >>> assert len(intensities) == len(chromatogram.get_time_array())
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() )
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0)

