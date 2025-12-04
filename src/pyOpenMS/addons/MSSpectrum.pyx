

    def get_mz_array(MSSpectrum self):
        """
        Get the m/z values of the spectrum as a numpy array.

        Returns:
            np.ndarray: A 1D numpy array (float64) containing the m/z values
                       for each peak in the spectrum.

        Example:
            >>> spectrum = MSSpectrum()
            >>> mz_values = spectrum.get_mz_array()
            >>> print(f"m/z range: {mz_values.min():.2f} - {mz_values.max():.2f}")
        """
        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef size_t n = spec_.size()

        if n == 0:
            return np.empty(0, dtype=np.float64)

        cdef np.ndarray[np.float64_t, ndim=1] mzs = np.empty(n, dtype=np.float64)
        cdef size_t i
        for i in range(n):
            mzs[i] = deref(spec_)[i].getMZ()

        return mzs

    def get_intensity_array(MSSpectrum self):
        """
        Get the intensity values of the spectrum as a numpy array.

        Returns:
            np.ndarray: A 1D numpy array (float32) containing the intensity values
                       for each peak in the spectrum.

        Example:
            >>> spectrum = MSSpectrum()
            >>> intensities = spectrum.get_intensity_array()
            >>> print(f"Total ion current: {intensities.sum():.2f}")
        """
        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef size_t n = spec_.size()

        if n == 0:
            return np.empty(0, dtype=np.float32)

        cdef np.ndarray[np.float32_t, ndim=1] intensities = np.empty(n, dtype=np.float32)
        cdef size_t i
        for i in range(n):
            intensities[i] = deref(spec_)[i].getIntensity()

        return intensities

    def get_peaks(self):
        """Cython signature: numpy_vector, numpy_vector get_peaks()

        Will return a tuple of two numpy arrays (m/z, intensity) corresponding
        to the peaks in the MSSpectrum. Provides fast access to peaks.

        Returns:
            tuple: A tuple of (mz_array, intensity_array) where:
                - mz_array is np.ndarray[float64] of m/z values
                - intensity_array is np.ndarray[float32] of intensity values

        Example:
            >>> spectrum = MSSpectrum()
            >>> spectrum.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))
            >>> mz, intensities = spectrum.get_peaks()
            >>> print(f"Base peak m/z: {mz[intensities.argmax()]}")
        """
        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef size_t n = spec_.size()

        if n == 0:
            return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float32)

        cdef np.ndarray[np.float64_t, ndim=1] mzs = np.empty(n, dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities = np.empty(n, dtype=np.float32)

        # Optimized: use direct indexing instead of iterator
        cdef size_t i
        for i in range(n):
            mzs[i] = deref(spec_)[i].getMZ()
            intensities[i] = deref(spec_)[i].getIntensity()

        return mzs, intensities

    def set_peaks(self, peaks):
        """Cython signature: set_peaks((numpy_vector, numpy_vector))
        
        Takes a tuple or list of two arrays (m/z, intensity) and populates the
        MSSpectrum. The arrays can be numpy arrays (faster).
        """

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"

        mzs, intensities = peaks
        assert len(mzs) == len(intensities), "Input vectors for set_peaks need to have the same length (mz and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float32 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(mzs, intensities)
        elif isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float64 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(mzs, intensities)
        else:
            self._set_peaks_orig(mzs, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef double intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_orig(self, mzs, intensities):


        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(mzs)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(mzs)

        for i in range(N):
            mz = mzs[i]
            intensity = intensities[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()

    def intensityInRange(self, float mzmin, float mzmax):

        cdef double I

        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef int N = spec_.size()

        I = 0.0
        for i in range(N):
                if deref(spec_)[i].getMZ() >= mzmin:
                    break

        cdef _Peak1D * p
        for j in range(i, N):
                p = address(deref(spec_)[j])
                if p.getMZ() > mzmax:
                    break
                I += p.getIntensity()

        return I

    def getIMData(self):
        """
        Get the position of ion mobility data array and its unit.

        Returns:
            tuple: (index, unit) where index is the position in FloatDataArrays
                   and unit is the DriftTimeUnit enum value.

        Raises:
            Exception: If no ion mobility data is present. Use containsIMData() first.

        Example:
            >>> if spectrum.containsIMData():
            ...     idx, unit = spectrum.getIMData()
            ...     im_array = spectrum.getFloatDataArrays()[idx]
        """
        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()

        pos = r.first
        unit = <int>r.second

        return (pos, unit)

    def get_drift_time_array(self):
        """
        Get the ion mobility drift time array as a numpy array (copy).

        This is a convenience method that retrieves the ion mobility data
        from the FloatDataArrays and returns it as a numpy array.

        Returns:
            np.ndarray or None: A 1D numpy array (float32) containing drift time
                               values for each peak, or None if no IM data present.

        Example:
            >>> spectrum = MSSpectrum()
            >>> drift_times = spectrum.get_drift_time_array()
            >>> if drift_times is not None:
            ...     print(f"Drift time range: {drift_times.min():.2f} - {drift_times.max():.2f}")
        """
        if not self.containsIMData():
            return None

        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()
        cdef size_t pos = r.first

        cdef libcpp_vector[_FloatDataArray] fdas = self.inst.get().getFloatDataArrays()
        cdef _FloatDataArray * fda = &fdas[pos]
        cdef size_t n = fda.size()

        if n == 0:
            return np.empty(0, dtype=np.float32)

        cdef np.ndarray[np.float32_t, ndim=1] result = np.empty(n, dtype=np.float32)
        cdef size_t i
        for i in range(n):
            result[i] = deref(fda)[i]

        return result

    def get_drift_time_array_mv(self):
        """
        Get the ion mobility drift time array as a memory view (no copy).

        This method provides direct access to the underlying drift time data
        without copying, which is more memory efficient for large datasets.

        Returns:
            memoryview or None: A memory view of drift time values, or None if
                               no IM data is present or array is empty.

        Warning:
            The returned memory view refers directly to the underlying data in
            a FloatDataArray. You must keep a reference to the FloatDataArray
            (via getFloatDataArrays()) to ensure the data remains valid.

            For safer access, use get_drift_time_array() which returns a copy.

        Example:
            >>> if spectrum.containsIMData():
            ...     # Keep reference to data arrays to prevent garbage collection
            ...     fdas = spectrum.getFloatDataArrays()
            ...     idx, unit = spectrum.getIMData()
            ...     drift_mv = spectrum.get_drift_time_array_mv()
            ...     total = sum(drift_mv)
        """
        if not self.containsIMData():
            return None

        # Get the position of IM data
        pos, unit = self.getIMData()
        fdas = self.getFloatDataArrays()

        if pos >= len(fdas):
            return None

        # Use FloatDataArray's get_data_mv() which returns a memory view
        return fdas[pos].get_data_mv()

    def get_drift_time_unit(self):
        """
        Get the drift time unit for ion mobility data.

        Returns:
            int or None: The DriftTimeUnit enum value, or None if no IM data present.
                        Values: 0=NONE, 1=MILLISECOND, 2=VSSC, 3=FAIMS_COMPENSATION_VOLTAGE

        Example:
            >>> unit = spectrum.get_drift_time_unit()
            >>> if unit == 1:  # DriftTimeUnit.MILLISECOND
            ...     print("Drift time is in milliseconds")
        """
        if not self.containsIMData():
            return None

        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()
        return <int>r.second

    def __len__(self):
        """Return the number of peaks in the spectrum."""
        return self.inst.get().size()
