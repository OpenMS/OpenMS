cimport numpy as np
import numpy as np






    def getMetaValues(self):
        """Cython signature: dict getMetaValues()

        Returns all meta values in a Python dictionary.
        """
        # cdef _DataValue _c_value
        # cdef DataValue _value 
        # cdef int _type

        mmap = {}
        cdef libcpp_vector[_String] keys 
        cdef object py_result

        # get keys for iteration
        self.inst.get().getKeys(keys)

        cdef libcpp_vector[_String].iterator k_it = keys.begin()
        while k_it != keys.end():
            # easy approach: call Python fxn
            py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            py_result = self.getMetaValue(py_str)
            # # hard approach: do it ourselves
            # _c_value = self.inst.get().getMetaValue(deref(k_it))
            # py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            # _type = _c_value.valueType()
            # _value = DataValue.__new__(DataValue)
            # _value.inst = shared_ptr[_DataValue](new _DataValue(_c_value))
            # if _type == DataType.STRING_VALUE:
            #     py_result = _value.toString()
            # elif _type == DataType.INT_VALUE:
            #     py_result = _value.toInt()
            # elif _type == DataType.DOUBLE_VALUE:
            #     py_result = _value.toDouble()
            # elif _type == DataType.INT_LIST:
            #     py_result = _value.toIntList()
            # elif _type == DataType.DOUBLE_LIST:
            #     py_result = _value.toDoubleList()
            # elif _type == DataType.STRING_LIST:
            #     py_result = _value.toStringList()
            # elif _type == DataType.EMPTY_VALUE:
            #     py_result = None
            # else:
            #     raise Exception("DataValue instance has invalid value type %d" % _type)

            mmap[ py_str ] = py_result
            inc(k_it)

        return mmap

    def setMetaValues(self, dict mmap):
        """Cython signature: setMetaValues(dict values)

        Sets the meta values given in the Python dictionary.
        """
        self.inst.get().clearMetaInfo() # ensure its empty first
        for k, v in mmap.iteritems():
            # self.inst.get().setMetaValue(deref((convString(k)).get()), deref(DataValue(v).inst.get()))
            self.setMetaValue(k, v)








    def get_peaks(self):
        """Cython signature: numpy_vector, numpy_vector get_peaks()
        
        Will return a tuple of two numpy arrays (m/z, intensity) corresponding
        to the peaks in the MSSpectrum. Provides fast access to peaks.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mzs
        mzs = np.zeros( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.zeros( (n,), dtype=np.float32)
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

        spec_.clear(0) # empty vector , keep meta data
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

        spec_.clear(0) # empty vector , keep meta data
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

        spec_.clear(0) # empty vector , keep meta data
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

        cdef int n
        cdef double I

        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef int N = spec_.size()

        I = 0
        for i in range(N):
                if deref(spec_)[i].getMZ() >= mzmin:
                    break

        cdef _Peak1D * p
        for j in range(i, N):
                p = address(deref(spec_)[i])
                if p.getMZ() > mzmax:
                    break
                I += p.getIntensity()

        return I

