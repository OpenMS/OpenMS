

    def get_mz_array(self):
        """
        get_mz_array(self: OSSpectrum) -> np.ndarray

        Get the m/z array of the spectrum as a numpy array (copy).

        Returns:
            np.ndarray: A 1D numpy array (float64) containing the mass-to-charge ratio values
                       for each peak in the spectrum.

        Example:
            >>> spectrum = OSSpectrum()
            >>> mz_values = spectrum.get_mz_array()
            >>> print(f"m/z range: {min(mz_values)} - {max(mz_values)}")
            >>> # Safe to modify the returned array
            >>> mz_values += 0.001  # This won't affect the original data
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_mz_array_mv(self):
        """
        get_mz_array_mv(self: OSSpectrum) -> Optional[memoryview]

        Get the m/z array of the spectrum as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A memory view of the mass-to-charge ratio values
                       for each peak in the spectrum, or None if empty.

        Warning:
            DO NOT store the returned memory view for later use after the
            spectrum object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> spectrum = OSSpectrum()
            >>> mz_mv = spectrum.get_mz_array_mv()
            >>> # CORRECT: Use immediately
            >>> min_mz = min(mz_mv)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_mv = spectrum.get_mz_array_mv()
            >>> # del spectrum  # Memory view now points to invalid memory!
            >>> # print(stored_mv[0])  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr

    def get_intensity_array(self):
        """
        get_intensity_array(self: OSSpectrum) -> np.ndarray

        Get the intensity array of the spectrum as a numpy array (copy).

        Returns:
            np.ndarray: A 1D numpy array (float64) containing the intensity values
                       for each peak in the spectrum.

        Example:
            >>> spectrum = OSSpectrum()
            >>> intensities = spectrum.get_intensity_array()
            >>> print(f"Base peak intensity: {max(intensities)}")
            >>> # Safe to modify the returned array
            >>> normalized = intensities / max(intensities)  # This won't affect the original data
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_intensity_array_mv(self):
        """
        get_intensity_array_mv(self: OSSpectrum) -> Optional[memoryview]

        Get the intensity array of the spectrum as a memory view (no copy).

        This method provides direct access to the underlying data without copying,
        which is more memory efficient for large datasets.

        Returns:
            memoryview: A memory view of the intensity values
                       for each peak in the spectrum, or None if empty.

        Warning:
            DO NOT store the returned memory view for later use after the
            spectrum object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> spectrum = OSSpectrum()
            >>> intensities_mv = spectrum.get_intensity_array_mv()
            >>> # CORRECT: Use immediately for calculations
            >>> total_intensity = sum(intensities_mv)
            >>>
            >>> # WRONG: Don't do this!
            >>> # stored_mv = spectrum.get_intensity_array_mv()
            >>> # del spectrum  # Memory view now points to invalid memory!
            >>> # print(sum(stored_mv))  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr

    def get_drift_time_array(self):
        """
        get_drift_time_array(self: OSSpectrum) -> Optional[np.ndarray]

        Get the drift time array of the spectrum as a numpy array (copy).

        This method is used for ion mobility spectrometry data where each peak
        has an associated drift time value.

        Returns:
            np.ndarray or None: A 1D numpy array (float64) containing the drift time values
                               for each peak in the spectrum, or None if no drift
                               time data is available.

        Example:
            >>> spectrum = OSSpectrum()
            >>> drift_times = spectrum.get_drift_time_array()
            >>> if drift_times is not None:
            ...     print(f"Drift time range: {min(drift_times)} - {max(drift_times)} ms")
            ... else:
            ...     print("No drift time data available")
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL:
            return None
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return np.empty(0, dtype=np.float64)
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return np.asarray(arr.copy())

    def get_drift_time_array_mv(self):
        """
        get_drift_time_array_mv(self: OSSpectrum) -> Optional[memoryview]

        Get the drift time array of the spectrum as a memory view (no copy).

        This method provides direct access to the underlying drift time data without copying,
        which is more memory efficient for large ion mobility datasets.

        Returns:
            memoryview or None: A memory view of the drift time values
                               for each peak in the spectrum, or None if no drift
                               time data is available.

        Warning:
            DO NOT store the returned memory view for later use after the
            spectrum object goes out of scope, as this will lead to
            undefined behavior and potential crashes.

        Example:
            >>> spectrum = OSSpectrum()
            >>> drift_mv = spectrum.get_drift_time_array_mv()
            >>> if drift_mv is not None:
            ...     # CORRECT: Use immediately
            ...     avg_drift = sum(drift_mv) / len(drift_mv)
            ...
            >>> # WRONG: Don't do this!
            >>> # stored_mv = spectrum.get_drift_time_array_mv()
            >>> # del spectrum  # Memory view now points to invalid memory!
            >>> # if stored_mv: print(stored_mv[0])  # Undefined behavior/crash
        """
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getDriftTimeArray()
        if _r.get() == NULL:
            return None
        cdef size_t n = _r.get().data.size()
        if n == 0:
            return None
        cdef double[::1] arr = <double [:n]>_r.get().data.data()
        return arr

    def get_data_arrays(self):
        """
        get_data_arrays(self: OSSpectrum) -> List[OSBinaryDataArray]

        Get all data arrays associated with the spectrum.

        This includes the m/z array, intensity array, drift time array (if present),
        and any additional meta data arrays that may be present in the spectrum.

        Returns:
            list: A list of OSBinaryDataArray objects containing all data arrays
                 associated with this spectrum.

        Example:
            >>> spectrum = OSSpectrum()
            >>> data_arrays = spectrum.get_data_arrays()
            >>> print(f"Number of data arrays: {len(data_arrays)}")
            >>> for i, array in enumerate(data_arrays):
            ...     data = array.get_data()
            ...     if data is not None:
            ...         print(f"Array {i}: {len(data)} data points")
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
        set_data_arrays(self: OSSpectrum, inp: List[OSBinaryDataArray]) -> None

        Set all data arrays for the spectrum.

        This method replaces all existing data arrays with the provided ones.
        The input should include m/z array, intensity array, and any additional
        meta data arrays (such as drift time arrays for ion mobility data).

        Args:
            inp (list): A list of OSBinaryDataArray objects to set as the
                       spectrum's data arrays.

        Raises:
            AssertionError: If inp is not a list or contains elements that are
                           not OSBinaryDataArray instances.

        Example:
            >>> spectrum = OSSpectrum()
            >>> mz_array = OSBinaryDataArray()
            >>> intensity_array = OSBinaryDataArray()
            >>> # ... populate arrays with data ...
            >>> spectrum.set_data_arrays([mz_array, intensity_array])
        """
        assert isinstance(inp, list) and all(isinstance(ele_inp, OSBinaryDataArray) for ele_inp in inp), 'Input has to be a list of elements of type OSBinaryDataArray'


        cdef libcpp_vector[ shared_ptr[_OSBinaryDataArray] ]  _r = self.inst.get().getDataArrays()
        _r.clear()

        cdef OSBinaryDataArray rv
        for rv in inp:
            _r.push_back(rv.inst)

        self.inst.get().setDataArrays(_r)

    def set_mz_array(self, list data):
        """
        set_mz_array(self: OSSpectrum, data: List[float]) -> None

        Set the m/z array (mass-to-charge ratios) for the spectrum.

        Args:
            data (list): A list of numeric values representing the m/z values
                        for each peak in the spectrum. Values should typically be
                        in ascending order for proper spectrum processing.

        Raises:
            AssertionError: If data is not a list.

        Example:
            >>> spectrum = OSSpectrum()
            >>> mz_values = [100.05, 200.10, 300.15, 400.20, 500.25]
            >>> spectrum.set_mz_array(mz_values)
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() )
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setMZArray(v0)


    def set_intensity_array(self, list data):
        """
        set_intensity_array(self: OSSpectrum, data: List[float]) -> None

        Set the intensity array for the spectrum.

        Args:
            data (list): A list of numeric values representing the intensity values
                        for each peak in the spectrum. The length should match the
                        m/z array length.

        Raises:
            AssertionError: If data is not a list.

        Example:
            >>> spectrum = OSSpectrum()
            >>> intensities = [1000.0, 5000.0, 3000.0, 2000.0, 1500.0]
            >>> spectrum.set_intensity_array(intensities)
            >>> # Make sure length matches m/z array
            >>> assert len(intensities) == len(spectrum.get_mz_array())
        """
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() )
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0)
