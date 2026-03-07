cimport numpy as np
import numpy as np




    def get_peaks(self):

        cdef _Mobilogram * mobilogram_ = self.inst.get()

        cdef unsigned int n = mobilogram_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mobilities
        mobilities = np.zeros( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.zeros( (n,), dtype=np.float32)
        cdef _MobilityPeak1D p

        cdef libcpp_vector[_MobilityPeak1D].iterator it = mobilogram_.begin()
        cdef int i = 0
        while it != mobilogram_.end():
            mobilities[i] = deref(it).getMobility()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return mobilities, intensities

    def set_peaks(self, peaks):

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mobility and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mobility and intensity vector)"

        mobilities, intensities = peaks
        assert len(mobilities) == len(intensities), "Input vectors for set_peaks need to have the same length (mobility and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(mobilities, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mobilities.dtype == np.float64 and intensities.dtype == np.float32 and \
          mobilities.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(mobilities, intensities)
        elif isinstance(mobilities, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mobilities.dtype == np.float64 and intensities.dtype == np.float64 and \
          mobilities.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(mobilities, intensities)
        else:
            self._set_peaks_orig(mobilities, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_mb not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _Mobilogram * mobilogram_ = self.inst.get()

        mobilogram_.resize(0) # empty vector, keep meta data and data arrays
        mobilogram_.reserve(<int>len(data_mb)) # allocate space for incoming data
        cdef _MobilityPeak1D p = _MobilityPeak1D()
        cdef double mb
        cdef double intensity
        cdef int N
        N = len(data_mb)

        for i in range(N):
            mb = data_mb[i]
            intensity = data_i[i]
            p.setMobility(<double>mb)
            p.setIntensity(<float>intensity)
            mobilogram_.push_back(p)

        mobilogram_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_mb not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _Mobilogram * mobilogram_ = self.inst.get()

        mobilogram_.resize(0) # empty vector, keep meta data and data arrays
        mobilogram_.reserve(<int>len(data_mb)) # allocate space for incoming data
        cdef _MobilityPeak1D p = _MobilityPeak1D()
        cdef double mb
        cdef float intensity
        cdef int N
        N = len(data_mb)

        for i in range(N):
            mb = data_mb[i]
            intensity = data_i[i]
            p.setMobility(<double>mb)
            p.setIntensity(<float>intensity)
            mobilogram_.push_back(p)

        mobilogram_.updateRanges()


    def _set_peaks_orig(self, mobilities, intensities):


        cdef _Mobilogram * mobilogram_ = self.inst.get()

        mobilogram_.resize(0) # empty vector, keep meta data and data arrays
        mobilogram_.reserve(<int>len(mobilities)) # allocate space for incoming data
        cdef _MobilityPeak1D p = _MobilityPeak1D()
        cdef double mb
        cdef float intensity
        cdef int N
        N = len(mobilities)

        for i in range(N):
            mb = mobilities[i]
            intensity = intensities[i]
            p.setMobility(<double>mb)
            p.setIntensity(<float>intensity)
            mobilogram_.push_back(p)

        mobilogram_.updateRanges()

    def __len__(self):
        """Return the number of peaks in the mobilogram."""
        return self.inst.get().size()

    def __str__(self):
        """
        Return a string representation of the Mobilogram object.
        Delegates to __repr__ for consistency.
        """
        return self.__repr__()

    def __repr__(self):
        """
        Return a string representation of the Mobilogram object.

        Returns key properties in a readable format:
        Mobilogram(num_peaks=100, rt=123.45, mobility_range=[0.5, 1.5])
        """
        cdef size_t num_peaks = self.inst.get().size()

        parts = []
        parts.append(f"num_peaks={num_peaks}")

        # Add RT
        rt = self.getRT()
        if rt > 0:
            parts.append(f"rt={rt:.2f}")

        # Add drift time unit
        unit_str = self.getDriftTimeUnitAsString()
        if unit_str:
            unit_decoded = unit_str.decode('utf-8') if isinstance(unit_str, bytes) else str(unit_str)
            if unit_decoded and unit_decoded != 'none':
                parts.append(f"unit='{unit_decoded}'")

        # Add mobility range if there are peaks
        if num_peaks > 0:
            mobilities, _ = self.get_peaks()
            parts.append(f"mobility_range=[{mobilities[0]:.4f}, {mobilities[-1]:.4f}]")

        return f"Mobilogram({', '.join(parts)})"

    def get_df_columns(self, columns='default'):
        """Returns a list of column names that get_df() would produce for this mobilogram.

        Useful for discovering available columns before export.

        Args:
            columns (str): 'default' for standard columns, 'all' for all available
                          columns including non-default ones.

        Returns:
            list: List of column name strings.

        Example:
            >>> cols = mobilogram.get_df_columns()
            ['mobility', 'intensity', 'rt', 'drift_time_unit']
        """
        # Default columns (Mobilogram doesn't support meta values)
        return ['mobility', 'intensity', 'rt', 'drift_time_unit']

    def get_data_dict(self, columns=None):
        """Returns a dictionary of NumPy arrays with mobility, intensities, and metadata.

        This method extracts mobilogram data including peaks
        into a dictionary format suitable for conversion to a pandas DataFrame.

        Args:
            columns (list or None): List of column names to include. If None, includes
                                   all default columns. Use get_df_columns() to see
                                   all available columns.

        Returns:
            dict: Dictionary with requested columns as keys and numpy arrays as values.
                - 'mobility': numpy array of mobility values (float64)
                - 'intensity': numpy array of intensity values (float32)
                - 'rt': retention time (float64)
                - 'drift_time_unit': drift time unit string

        Example:
            >>> data = mobilogram.get_data_dict()
            >>> data = mobilogram.get_data_dict(columns=['mobility', 'intensity'])
        """
        # Get peak data
        cdef np.ndarray[np.float64_t, ndim=1] mobilities
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        mobilities, intensities = self.get_peaks()
        cnt = len(mobilities)

        # Determine which columns to include
        if columns is not None:
            requested = set(columns)
        else:
            requested = None  # None means include all defaults

        def want(col):
            """Check if a default column should be included."""
            return requested is None or col in requested

        data_dict = {}

        # Core peak data
        if want('mobility'):
            data_dict['mobility'] = mobilities
        if want('intensity'):
            data_dict['intensity'] = intensities

        # RT
        if want('rt'):
            data_dict['rt'] = np.full(cnt, self.getRT(), dtype=np.float64)

        # Drift time unit
        if want('drift_time_unit'):
            unit_str = self.getDriftTimeUnitAsString()
            unit_decoded = unit_str.decode('utf-8') if isinstance(unit_str, bytes) else str(unit_str)
            data_dict['drift_time_unit'] = np.full(cnt, unit_decoded, dtype='U50')

        return data_dict
