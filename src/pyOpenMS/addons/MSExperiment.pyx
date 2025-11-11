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
