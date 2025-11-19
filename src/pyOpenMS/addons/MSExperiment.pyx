import numpy as np


    def get2DPeakData(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakData(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""        

        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        exp_.get2DPeakDataPerSpectrum(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty)



    def get2DPeakDataLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakDataLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        exp_.get2DPeakData(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap))
    
    def get2DPeakDataIM(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIM(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        cdef libcpp_vector[libcpp_vector[float]] ion_mobility
        exp_.get2DPeakDataIMPerSpectrum(min_rt, max_rt, min_mz, max_mz,  ms_level, rt, mz, inty, ion_mobility)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_ion = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap
        cdef ArrayWrapperFloat ion_mobility_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            ion_mobility_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            ion_mobility_wrap.set_data(ion_mobility[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)
            all_ion[i] = np.frombuffer(ion_mobility_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty, all_ion)

    def get2DPeakDataIMLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIMLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        cdef libcpp_vector[float] ion_mobility
        exp_.get2DPeakDataIM(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty, ion_mobility)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat ion_mobility_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)
        ion_mobility_wrap.set_data(ion_mobility)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap), np.asarray(ion_mobility_wrap))

    def getMSLevels(self):
        """Cython signature: list[int] getMSLevels()"""
        cdef libcpp_vector[unsigned int] _r = self.inst.get().getMSLevels()
        cdef libcpp_vector[unsigned int].iterator it__r = _r.begin()
        cdef list result = []
        while it__r != _r.end():
            result.append(deref(it__r))
            inc(it__r)
        return result

    def getChromatogram(self, id_):
        """Cython signature: `MSChromatogram getChromatogram(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrChromatograms(), 'Requested chromatogram %s does not exist, there are only %s chromatograms' % (id_, self.getNrChromatograms() )
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>id_)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result

    def getSpectrum(self, id_):
        """Cython signature: `MSSpectrum getSpectrum(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrSpectra(), 'Requested spectrum %s does not exist, there are only %s spectra' % (id_, self.getNrSpectra() )
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>id_)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result

    @property
    def ms_level(self):
        """
        Extract MS levels from all spectra in the experiment.

        Returns a NumPy array containing the MS level for each spectrum in the
        MSExperiment. The array indices correspond to spectrum indices, allowing
        for direct mapping between spectra and their MS levels.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of unsigned integers (uint32) containing the MS level
            of each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> ms_levels = exp.ms_level
        >>> print(ms_levels)
        [1 2 2 1 2 2 ...]
        >>> # Find all MS2 spectra indices
        >>> ms2_indices = np.where(ms_levels == 2)[0]

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - MS levels are typically 1 (MS1) or 2 (MS2), but higher levels are possible
        - The returned array is a standard NumPy array suitable for vectorized operations
        - This method is more efficient than calling `getSpectrum(i).getMSLevel()`
          repeatedly for large datasets

        See Also
        --------
        getMSLevels : Get unique MS levels present in the experiment
        getSpectrum : Retrieve individual spectrum by index
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[unsigned int] ms_levels
        ms_levels.reserve(n_spectra)
        
        cdef size_t i
        for i in range(n_spectra):
            ms_levels.push_back(self.inst.get().getSpectrum(i).getMSLevel())
        
        return np.asarray(ms_levels)

    @property
    def rt(self):
        """
        Access retention times with support for indexed assignment.
        
        Returns an RTAccessor object that behaves like a numpy array but
        supports indexed assignment operations. This enables direct updates
        to retention times using boolean masks, integer arrays, or slices.

        Returns
        -------
        RTAccessor
            An accessor object that supports numpy-style operations and
            indexed assignment. Can be converted to a regular numpy array
            using np.array(exp.rt).

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>>
        >>> # Read operations (backward compatible)
        >>> rts = exp.rt  # Returns RTAccessor
        >>> print(rts)  # Works like numpy array
        [0.5 1.2 1.9 2.6 ...]
        >>> rts_array = np.array(exp.rt)  # Convert to numpy array
        >>>
        >>> # Indexed assignment (new feature)
        >>> exp.rt[exp.ms_level == 2] = [120.5, 180.3]  # Boolean indexing
        >>> exp.rt[[0, 2, 4]] = [100.0, 200.0, 300.0]  # Integer array indexing
        >>> exp.rt[1:4] = [150.0, 175.0, 200.0]  # Slice indexing
        >>> exp.rt[5] = 250.0  # Single value assignment
        >>>
        >>> # Works with filtering
        >>> mask = (exp.rt >= 1.0) & (exp.rt <= 2.0)
        >>> filtered_exp = exp[mask]

        Notes
        -----
        - RTAccessor implements the numpy array protocol for full compatibility
        - All existing code using exp.rt should continue to work without changes
        - Retention times are in seconds
        - For performance-critical code, convert to numpy array once and reuse
        - The rt.setter is still available for setting all values at once

        See Also
        --------
        getSpectrum : Retrieve individual spectrum by index
        rt.setter : Set all retention times at once
        """
        return RTAccessor(self)

    # @rt.setter
    # def rt(self, values):
    #     """
    #     Set retention times for all spectra.
        
    #     Parameters
    #     ----------
    #     values : array-like
    #         Array of retention times (in seconds) with length equal to number of spectra
            
    #     Raises
    #     ------
    #     ValueError
    #         If length of values doesn't match number of spectra
    #     """
    #     values_np = np.asarray(values, dtype=np.float64)
    #     cdef size_t n_spectra = self.inst.get().getNrSpectra()
        
    #     if len(values_np) != n_spectra:
    #         raise ValueError(f"Length of values ({len(values_np)}) must match number of spectra ({n_spectra})")
        
    #     # Get copy of all spectra
    #     cdef libcpp_vector[_MSSpectrum] spectra = self.inst.get().getSpectra()
        
    #     # Modify the copy
    #     cdef size_t i
    #     for i in range(n_spectra):
    #         spectra[i].setRT(values_np[i])
        
    #     # Set modified spectra back
    #     self.inst.get().setSpectra(spectra)

    @property
    def tic(self):
        """
        Extract total ion current (TIC) from all spectra in the experiment.

        Calculates and returns a NumPy array containing the total ion current
        (sum of all peak intensities) for each spectrum in the MSExperiment.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the TIC of each spectrum.
            The length of the array equals the number of spectra in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> tics = exp.tic
        >>> print(tics)
        [1234.5 2345.6 3456.7 ...]
        >>> # Find spectra with high TIC
        >>> high_tic_mask = tics > np.percentile(tics, 75)
        >>> high_tic_spectra = exp[high_tic_mask]

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - TIC is calculated as the sum of all peak intensities in a spectrum
        - Useful for quality control and identifying high-signal spectra
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        max_intensity : Get maximum intensity values
        n_peaks : Get number of peaks per spectrum
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] tics
        tics.reserve(n_spectra)
        
        cdef size_t i
        for i in range(n_spectra):
            tics.push_back(self.inst.get().getSpectrum(i).calculateTIC())
        
        return np.asarray(tics)

    @property
    def n_peaks(self):
        """
        Extract number of peaks from all spectra in the experiment.

        Returns a NumPy array containing the peak count for each spectrum
        in the MSExperiment.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of unsigned integers containing the number of peaks
            in each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> peak_counts = exp.n_peaks
        >>> print(peak_counts)
        [123 234 345 ...]
        >>> # Find spectra with at least 100 peaks
        >>> rich_spectra = exp[peak_counts >= 100]
        >>> # Get average peaks per spectrum
        >>> avg_peaks = np.mean(peak_counts)

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Peak count includes all peaks in a spectrum, regardless of intensity
        - Empty spectra will have a count of 0
        - Useful for quality control and filtering

        See Also
        --------
        tic : Get total ion current values
        max_intensity : Get maximum intensity values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[size_t] peak_counts
        peak_counts.reserve(n_spectra)
        
        cdef size_t i
        for i in range(n_spectra):
            peak_counts.push_back(self.inst.get().getSpectrum(i).size())
        
        return np.asarray(peak_counts)

    @property
    def max_intensity(self):
        """
        Extract maximum intensity from all spectra in the experiment.

        Returns a NumPy array containing the maximum peak intensity for each
        spectrum in the MSExperiment. For empty spectra, returns 0.0.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the maximum intensity
            of each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> max_ints = exp.max_intensity
        >>> print(max_ints)
        [9876.5 8765.4 7654.3 ...]
        >>> # Find spectra with high max intensity
        >>> intense_spectra = exp[max_ints > 5000]
        >>> # Normalize by max intensity
        >>> normalized = max_ints / np.max(max_ints)

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Empty spectra return 0.0 as maximum intensity
        - This is computed by finding the peak with highest intensity in each spectrum
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        min_intensity : Get minimum intensity values
        tic : Get total ion current values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] max_intensities
        max_intensities.reserve(n_spectra)
        
        cdef size_t i, j
        cdef double max_int
        cdef size_t spec_size
        for i in range(n_spectra):
            spec_size = self.inst.get().getSpectrum(i).size()
            if spec_size == 0:
                max_intensities.push_back(0.0)
            else:
                max_int = self.inst.get().getSpectrum(i)[0].getIntensity()
                for j in range(1, spec_size):
                    if self.inst.get().getSpectrum(i)[j].getIntensity() > max_int:
                        max_int = self.inst.get().getSpectrum(i)[j].getIntensity()
                max_intensities.push_back(max_int)
        
        return np.asarray(max_intensities)

    @property
    def min_intensity(self):
        """
        Extract minimum intensity from all spectra in the experiment.

        Returns a NumPy array containing the minimum peak intensity for each
        spectrum in the MSExperiment. For empty spectra, returns 0.0.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the minimum intensity
            of each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> min_ints = exp.min_intensity
        >>> print(min_ints)
        [10.5 12.3 8.7 ...]
        >>> # Find spectra with noise floor above threshold
        >>> clean_spectra = exp[min_ints > 50]
        >>> # Get intensity range per spectrum
        >>> intensity_range = exp.max_intensity - exp.min_intensity

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Empty spectra return 0.0 as minimum intensity
        - This is computed by finding the peak with lowest intensity in each spectrum
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        max_intensity : Get maximum intensity values
        tic : Get total ion current values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] min_intensities
        min_intensities.reserve(n_spectra)
        
        cdef size_t i, j
        cdef double min_int
        cdef size_t spec_size
        for i in range(n_spectra):
            spec_size = self.inst.get().getSpectrum(i).size()
            if spec_size == 0:
                min_intensities.push_back(0.0)
            else:
                min_int = self.inst.get().getSpectrum(i)[0].getIntensity()
                for j in range(1, spec_size):
                    if self.inst.get().getSpectrum(i)[j].getIntensity() < min_int:
                        min_int = self.inst.get().getSpectrum(i)[j].getIntensity()
                min_intensities.push_back(min_int)
        
        return np.asarray(min_intensities)

    @property
    def mz_min(self):
        """
        Extract minimum m/z value from all spectra in the experiment.

        Returns a NumPy array containing the minimum m/z value for each
        spectrum in the MSExperiment. For empty spectra, returns 0.0.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the minimum m/z value
            of each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> min_mzs = exp.mz_min
        >>> print(min_mzs)
        [100.5 150.2 200.8 ...]
        >>> # Find spectra covering low m/z range
        >>> low_mz_spectra = exp[min_mzs < 200]
        >>> # Get m/z range per spectrum
        >>> mz_range = exp.mz_max - exp.mz_min

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Empty spectra return 0.0 as minimum m/z
        - This is computed by finding the peak with lowest m/z in each spectrum
        - Assumes spectrum is sorted by m/z (which is typical)
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        mz_max : Get maximum m/z values
        rt : Get retention time values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] min_mzs
        min_mzs.reserve(n_spectra)
        
        cdef size_t i
        for i in range(n_spectra):
            if self.inst.get().getSpectrum(i).size() == 0:
                min_mzs.push_back(0.0)
            else:
                min_mzs.push_back(self.inst.get().getSpectrum(i)[0].getMZ())
        
        return np.asarray(min_mzs)

    @property
    def mz_max(self):
        """
        Extract maximum m/z value from all spectra in the experiment.

        Returns a NumPy array containing the maximum m/z value for each
        spectrum in the MSExperiment. For empty spectra, returns 0.0.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the maximum m/z value
            of each spectrum. The length of the array equals the number of spectra
            in the experiment.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with spectra ...
        >>> max_mzs = exp.mz_max
        >>> print(max_mzs)
        [1999.5 1850.2 2000.8 ...]
        >>> # Find spectra covering high m/z range
        >>> high_mz_spectra = exp[max_mzs > 1500]
        >>> # Get m/z range per spectrum
        >>> mz_range = exp.mz_max - exp.mz_min

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Empty spectra return 0.0 as maximum m/z
        - This is computed by finding the peak with highest m/z in each spectrum
        - Assumes spectrum is sorted by m/z (which is typical)
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        mz_min : Get minimum m/z values
        rt : Get retention time values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] max_mzs
        max_mzs.reserve(n_spectra)
        
        cdef size_t i
        cdef size_t spec_size
        for i in range(n_spectra):
            spec_size = self.inst.get().getSpectrum(i).size()
            if spec_size == 0:
                max_mzs.push_back(0.0)
            else:
                max_mzs.push_back(self.inst.get().getSpectrum(i)[spec_size - 1].getMZ())
        
        return np.asarray(max_mzs)

    @property
    def drift_time(self):
        """
        Extract drift time from all spectra in the experiment.

        Returns a NumPy array containing the ion mobility drift time for each
        spectrum in the MSExperiment. For spectra without drift time information,
        returns NaN.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of doubles containing the drift time of each spectrum.
            The length of the array equals the number of spectra in the experiment.
            NaN values indicate missing drift time information.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with IM spectra ...
        >>> drift_times = exp.drift_time
        >>> print(drift_times)
        [10.5 11.2 nan 12.3 ...]
        >>> # Find spectra with drift time data
        >>> has_dt = ~np.isnan(drift_times)
        >>> im_spectra = exp[has_dt]
        >>> # Filter by drift time range
        >>> dt_mask = (drift_times >= 10) & (drift_times <= 15)
        >>> filtered = exp[dt_mask]

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - NaN indicates that drift time is not set for a spectrum
        - Drift time units can vary; check `drift_time_unit` property
        - Drift times may be stored directly as a spectrum attribute
        - The returned array is a standard NumPy array suitable for vectorized operations

        See Also
        --------
        drift_time_unit : Get drift time units
        rt : Get retention time values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[double] drift_times
        drift_times.reserve(n_spectra)
        
        cdef size_t i
        cdef double dt
        for i in range(n_spectra):
            dt = self.inst.get().getSpectrum(i).getDriftTime()
            # IMTypes::DRIFTTIME_NOT_SET is -1
            if dt < 0:
                drift_times.push_back(float('nan'))
            else:
                drift_times.push_back(dt)
        
        return np.asarray(drift_times)

    @drift_time.setter
    def drift_time(self, values):
        """
        Set drift times for all spectra.
        
        Parameters
        ----------
        values : array-like
            Array of drift times with length equal to number of spectra.
            Use NaN or negative values to indicate drift time is not set.
            
        Raises
        ------
        ValueError
            If length of values doesn't match number of spectra
        """
        values_np = np.asarray(values, dtype=np.float64)
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        
        if len(values_np) != n_spectra:
            raise ValueError(f"Length of values ({len(values_np)}) must match number of spectra ({n_spectra})")
        
        # Get copy of all spectra
        cdef libcpp_vector[_MSSpectrum] spectra = self.inst.get().getSpectra()
        
        # Modify the copy
        cdef size_t i
        cdef double dt
        for i in range(n_spectra):
            dt = values_np[i]
            # Convert NaN to DRIFTTIME_NOT_SET (-1)
            if np.isnan(dt):
                dt = -1.0
            spectra[i].setDriftTime(dt)
        
        # Set modified spectra back
        self.inst.get().setSpectra(spectra)

    @property
    def drift_time_unit(self):
        """
        Extract drift time units from all spectra in the experiment.

        Returns a NumPy array containing the drift time unit as a string for each
        spectrum in the MSExperiment. Common units include 'millisecond',
        'inverse_reduced_ion_mobility', 'volt-second_per_square_centimeter', etc.

        Returns
        -------
        numpy.ndarray
            A 1D numpy array of strings (dtype=object) containing the drift time unit
            of each spectrum. The length of the array equals the number of spectra
            in the experiment. Empty strings indicate no drift time unit is set.

        Examples
        --------
        >>> exp = MSExperiment()
        >>> # ... load or populate experiment with IM spectra ...
        >>> dt_units = exp.drift_time_unit
        >>> print(dt_units)
        ['millisecond' 'millisecond' '' 'millisecond' ...]
        >>> # Check consistency of units
        >>> unique_units = np.unique(dt_units[dt_units != ''])
        >>> print(unique_units)
        ['millisecond']

        Notes
        -----
        - Returns an empty array if the experiment contains no spectra
        - Empty string indicates no drift time unit is set (DriftTimeUnit::NONE)
        - Common units in ion mobility MS include:
          * 'millisecond' - time-based drift time
          * 'inverse_reduced_ion_mobility' - 1/K0 values
          * 'volt-second_per_square_centimeter' - collision cross section related
        - The returned array is a standard NumPy object array

        See Also
        --------
        drift_time : Get drift time values
        """
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef np.ndarray dt_units = np.empty(n_spectra, dtype=object)
        
        cdef size_t i
        cdef _String dt_unit_str
        for i in range(n_spectra):
            dt_unit_str = self.inst.get().getSpectrum(i).getDriftTimeUnitAsString()
            dt_units[i] = <bytes>dt_unit_str.c_str()
            if isinstance(dt_units[i], bytes):
                dt_units[i] = dt_units[i].decode('utf-8')
        
        return dt_units

    @ms_level.setter
    def ms_level(self, values):
        """
        Set MS levels for all spectra.
        
        Parameters
        ----------
        values : array-like of int
            Array of MS levels with length equal to number of spectra.
            Typically 1 for MS1, 2 for MS2, etc.
            
        Raises
        ------
        ValueError
            If length of values doesn't match number of spectra
        """
        values_np = np.asarray(values, dtype=np.uint32)
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        
        if len(values_np) != n_spectra:
            raise ValueError(f"Length of values ({len(values_np)}) must match number of spectra ({n_spectra})")
        
        # Get copy of all spectra
        cdef libcpp_vector[_MSSpectrum] spectra = self.inst.get().getSpectra()
        
        # Modify the copy
        cdef size_t i
        for i in range(n_spectra):
            spectra[i].setMSLevel(<unsigned int>values_np[i])
        
        # Set modified spectra back
        self.inst.get().setSpectra(spectra)

    def __getitem__(self, key):
        cdef size_t n_spectra = self.inst.get().getNrSpectra()
        cdef libcpp_vector[size_t] indices
        cdef size_t i
        
        # Convert key to numpy array for easier handling
        try:
            key_array = np.asarray(key)
        except:
            raise TypeError("Index must be array-like (list, numpy array, etc.)")
        
        # Determine if this is boolean indexing or fancy indexing
        if key_array.dtype == np.bool_:
            # Boolean indexing
            # Validate boolean mask length
            if len(key_array) != n_spectra:
                raise IndexError(
                    f"Boolean index has wrong length: {len(key_array)} instead of {n_spectra}"
                )
            # Convert boolean mask to indices
            indices_np = np.where(key_array)[0]
        else:
            # Fancy indexing with integers
            indices_np = key_array.astype(np.intp)
            # Validate indices are within bounds
            if len(indices_np) > 0:
                if np.any(indices_np < 0) or np.any(indices_np >= n_spectra):
                    raise IndexError(
                        f"Index out of bounds. Valid range is [0, {n_spectra - 1}]"
                    )
        
        # Convert numpy indices to C++ vector
        indices.reserve(len(indices_np))
        for i in range(len(indices_np)):
            indices.push_back(<size_t>indices_np[i])
        
        # Create new MSExperiment
        cdef _MSExperiment * new_exp = new _MSExperiment()
        
        # Copy selected spectra to new experiment
        cdef _MSSpectrum spec
        for i in range(indices.size()):
            spec = self.inst.get().getSpectrum(indices[i])
            new_exp.addSpectrum(spec)
        
        # Create Python wrapper for the new experiment
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](new_exp)
        
        return py_result

# RTAccessor class for indexed assignment support
cdef class RTAccessor:
    """
    Proxy class for retention time access with indexed assignment support.
    
    This class enables numpy-style indexed assignment to retention times while
    maintaining backward compatibility with existing code. It intercepts write
    operations through __setitem__ and delegates read operations to numpy arrays.
    
    Examples
    --------
    >>> exp.rt[exp.ms_level == 2] = [120.5, 180.3]  # Boolean indexing
    >>> exp.rt[[0, 2, 4]] = [100.0, 200.0, 300.0]  # Integer array indexing
    >>> exp.rt[1:4] = [150.0, 175.0, 200.0]  # Slice indexing
    >>> rts = np.array(exp.rt)  # Convert to regular numpy array
    """
    cdef MSExperiment parent_exp
    
    def __init__(self, MSExperiment parent):
        """Initialize RTAccessor with parent MSExperiment."""
        self.parent_exp = parent
    
    def __getitem__(self, key):
        """
        Get retention time values for the given key.
        
        Supports all numpy indexing types: integers, slices, boolean masks,
        and integer arrays. Delegates to numpy array for read operations.
        
        Parameters
        ----------
        key : int, slice, array-like
            Index, slice, boolean mask, or integer array
            
        Returns
        -------
        float or np.ndarray
            Retention time value(s) at the specified indices
        """
        # Get current RT values as numpy array
        rt_values = self._get_rt_array()
        # Delegate to numpy for indexing
        return rt_values[key]
    
    def __setitem__(self, key, values):
        """
        Set retention time values for the given key.
        
        Intercepts write operations to update the underlying MSExperiment.
        Handles boolean masks, integer arrays, and slices.
        
        Parameters
        ----------
        key : int, slice, array-like
            Index, slice, boolean mask, or integer array
        values : float or array-like
            New retention time value(s)
            
        Raises
        ------
        ValueError
            If length of values doesn't match number of selected indices
        IndexError
            If indices are out of bounds
        """
        # Get current RT values
        rt_values = self._get_rt_array()
        
        # Convert key to indices
        if isinstance(key, (int, np.integer)):
            # Single integer index
            indices = np.array([key], dtype=np.intp)
            values_array = np.array([values], dtype=np.float64)
        elif isinstance(key, slice):
            # Slice indexing
            indices = np.arange(len(rt_values))[key]
            values_array = np.asarray(values, dtype=np.float64)
        else:
            # Array-like key (boolean mask or integer array)
            key_array = np.asarray(key)
            if key_array.dtype == np.bool_:
                # Boolean mask - convert to indices
                if len(key_array) != len(rt_values):
                    raise IndexError(
                        f"Boolean mask length ({len(key_array)}) must match "
                        f"number of spectra ({len(rt_values)})"
                    )
                indices = np.where(key_array)[0]
            else:
                # Integer array
                indices = key_array.astype(np.intp)
            values_array = np.asarray(values, dtype=np.float64)
        
        # Validate indices
        if len(indices) > 0:
            if np.any(indices < 0) or np.any(indices >= len(rt_values)):
                raise IndexError(
                    f"Index out of bounds. Valid range is [0, {len(rt_values) - 1}]"
                )
        
        # Handle broadcasting of single value
        if values_array.ndim == 0 or (values_array.ndim == 1 and len(values_array) == 1):
            # Broadcast single value to all selected indices
            if values_array.ndim == 0:
                values_array = np.full(len(indices), values_array.item(), dtype=np.float64)
            else:
                values_array = np.full(len(indices), values_array[0], dtype=np.float64)
        elif len(values_array) != len(indices):
            raise ValueError(
                f"Length of values ({len(values_array)}) must match "
                f"number of selected indices ({len(indices)})"
            )
        
        # Update RT values at selected indices
        rt_values[indices] = values_array
        
        # Set back to experiment using existing setter
        self.parent_exp.rt = rt_values
    
    def __array__(self):
        """
        Enable numpy.array() conversion.
        
        Returns
        -------
        np.ndarray
            Numpy array of all retention time values
        """
        return self._get_rt_array()
    
    def __len__(self):
        """
        Return number of spectra.
        
        Returns
        -------
        int
            Number of spectra in the experiment
        """
        return self.parent_exp.getNrSpectra()
    
    def _get_rt_array(self):
        """
        Internal method to get RT values as numpy array.
        
        Returns
        -------
        np.ndarray
            Array of retention time values
        """
        cdef size_t n_spectra = self.parent_exp.inst.get().getNrSpectra()
        cdef libcpp_vector[double] rts
        rts.reserve(n_spectra)
        
        cdef size_t i
        for i in range(n_spectra):
            rts.push_back(self.parent_exp.inst.get().getSpectrum(i).getRT())
        
        return np.asarray(rts)
    
    # Arithmetic operators for numpy compatibility
    def __add__(self, other):
        return self._get_rt_array() + other
    
    def __radd__(self, other):
        return other + self._get_rt_array()
    
    def __sub__(self, other):
        return self._get_rt_array() - other
    
    def __rsub__(self, other):
        return other - self._get_rt_array()
    
    def __mul__(self, other):
        return self._get_rt_array() * other
    
    def __rmul__(self, other):
        return other * self._get_rt_array()
    
    def __truediv__(self, other):
        return self._get_rt_array() / other
    
    def __rtruediv__(self, other):
        return other / self._get_rt_array()
    
    # Comparison operators
    def __eq__(self, other):
        return self._get_rt_array() == other
    
    def __ne__(self, other):
        return self._get_rt_array() != other
    
    def __lt__(self, other):
        return self._get_rt_array() < other
    
    def __le__(self, other):
        return self._get_rt_array() <= other
    
    def __gt__(self, other):
        return self._get_rt_array() > other
    
    def __ge__(self, other):
        return self._get_rt_array() >= other
    
    def __iter__(self):
        """Iterate over retention time values."""
        return iter(self._get_rt_array())
    
    def __repr__(self):
        """String representation of RTAccessor."""
        return f"RTAccessor({self._get_rt_array()})"
    
    def __str__(self):
        """String conversion of RTAccessor."""
        return str(self._get_rt_array())
