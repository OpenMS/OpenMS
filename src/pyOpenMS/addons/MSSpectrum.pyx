import numpy as _np



    def get_data_dict(self, export_meta_values=True):
        """Returns a dictionary of NumPy arrays with m/z, intensities, and metadata.

        This method extracts spectrum data including peaks, retention time, MS level,
        ion mobility data (if present), precursor information, and optional meta values
        into a dictionary format suitable for conversion to a pandas DataFrame.

        Args:
            export_meta_values (bool): Whether to include meta values in the output.
                                       Defaults to True.

        Returns:
            dict: Dictionary with the following keys:
                - 'mz': numpy array of m/z values (float64)
                - 'intensity': numpy array of intensity values (float32)
                - 'rt': numpy array of retention time values (float64)
                - 'ms_level': numpy array of MS level values (uint16)
                - 'native_id': numpy array of native ID strings
                - 'ion_mobility': numpy array of ion mobility values (if present)
                - 'ion_mobility_unit': numpy array of ion mobility unit strings
                - 'precursor_mz': numpy array of precursor m/z values
                - 'precursor_charge': numpy array of precursor charge values
                - Additional keys for each meta value (if export_meta_values=True)
        """
        # Get peak data using existing optimized method
        cdef _np.ndarray[_np.float64_t, ndim=1] mzs
        cdef _np.ndarray[_np.float32_t, ndim=1] intensities
        mzs, intensities = self.get_peaks()

        cnt = len(mzs)
        data_dict = {
            'mz': mzs,
            'intensity': intensities,
            'rt': _np.full(cnt, self.getRT(), dtype=_np.float64),
            'ms_level': _np.full(cnt, self.getMSLevel(), dtype=_np.uint16),
            'native_id': _np.full(cnt, self.getNativeID(), dtype='U100')
        }

        # Ion mobility handling
        if self.containsIMData():
            im_index, drift_time_unit = self.getIMData()
            im_arrays = self.getFloatDataArrays()

            if im_index >= 0 and im_index < im_arrays.size():
                data_dict['ion_mobility'] = _np.array(
                    [im_arrays[im_index][i] for i in range(cnt)],
                    dtype=_np.float64
                )
            else:
                data_dict['ion_mobility'] = _np.full(cnt, _np.nan, dtype=_np.float64)
            data_dict['ion_mobility_unit'] = _np.full(
                cnt,
                self.getDriftTimeUnitAsString(),
                dtype='U50'
            )
        else:
            data_dict['ion_mobility'] = _np.full(cnt, _np.nan, dtype=_np.float64)
            data_dict['ion_mobility_unit'] = _np.full(cnt, '', dtype='U1')

        # Precursor Info
        precursors = self.getPrecursors()
        if precursors.size() > 0:
            precursor = precursors[0]
            data_dict['precursor_mz'] = _np.full(cnt, precursor.getMZ(), dtype=_np.float64)
            data_dict['precursor_charge'] = _np.full(cnt, precursor.getCharge(), dtype=_np.int16)
        else:
            data_dict['precursor_mz'] = _np.full(cnt, _np.nan, dtype=_np.float64)
            data_dict['precursor_charge'] = _np.full(cnt, 0, dtype=_np.int16)

        # Metadata Handling - use Python type introspection
        if export_meta_values:
            mvs = []
            self.getKeys(mvs)
            for k in mvs:
                if not self.metaValueExists(k):
                    continue
                # Get meta value - autowrap returns Python types
                v = self.getMetaValue(k)
                k_str = k.decode() if isinstance(k, bytes) else k

                try:
                    if isinstance(v, int):
                        data_dict[k_str] = _np.full(cnt, v, dtype=_np.int64)
                    elif isinstance(v, float):
                        data_dict[k_str] = _np.full(cnt, v, dtype=_np.float64)
                    elif isinstance(v, str):
                        data_dict[k_str] = _np.full(cnt, v, dtype=f"U{max(len(v), 1)}")
                    elif isinstance(v, bool):
                        data_dict[k_str] = _np.full(cnt, v, dtype=bool)
                    else:
                        # Fallback for other types
                        data_dict[k_str] = _np.full(cnt, str(v), dtype='object')
                except Exception:
                    data_dict[k_str] = _np.full(cnt, str(v), dtype='object')

        return data_dict



    def get_mz_array(MSSpectrum self):
        """Cython signature: numpy_vector get_mz_array()

        Will return a numpy array corresponding
        to the mz values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef _np.ndarray[_np.float64_t, ndim=1] mzs
        mzs = _np.empty( (n,), dtype=_np.float64)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            inc(it)
            i += 1

        return mzs


    def get_intensity_array(MSSpectrum self):
        """Cython signature: numpy_vector get_intensity_array()

        Will return a numpy array corresponding
        to the intensity values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef _np.ndarray[_np.float32_t, ndim=1] intensities
        intensities = _np.empty( (n,), dtype=_np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return intensities

    def get_peaks(self):
        """Cython signature: numpy_vector, numpy_vector get_peaks()

        Will return a tuple of two numpy arrays (m/z, intensity) corresponding
        to the peaks in the MSSpectrum. Provides fast access to peaks.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef _np.ndarray[_np.float64_t, ndim=1] mzs
        mzs = _np.empty( (n,), dtype=_np.float64)
        cdef _np.ndarray[_np.float32_t, ndim=1] intensities
        intensities = _np.empty( (n,), dtype=_np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

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
        if isinstance(mzs, _np.ndarray) and isinstance(intensities, _np.ndarray) and \
          mzs.dtype == _np.float64 and intensities.dtype == _np.float32 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(mzs, intensities)
        elif isinstance(mzs, _np.ndarray) and isinstance(intensities, _np.ndarray) and \
          mzs.dtype == _np.float64 and intensities.dtype == _np.float64 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(mzs, intensities)
        else:
            self._set_peaks_orig(mzs, intensities)



    def _set_peaks_fast_dd(self, _np.ndarray[double, ndim=1, mode="c"] data_mz not None, _np.ndarray[double, ndim=1, mode="c"] data_i not None):

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


    def _set_peaks_fast_df(self, _np.ndarray[double, ndim=1, mode="c"] data_mz not None, _np.ndarray[float, ndim=1, mode="c"] data_i not None):

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

        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()

        pos = r.first
        unit = <int>r.second

        return (pos, unit)

    def __len__(self):
        """Return the number of peaks in the spectrum."""
        return self.inst.get().size()
