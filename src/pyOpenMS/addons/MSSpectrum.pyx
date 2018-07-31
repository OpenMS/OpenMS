cimport numpy as np
import numpy as np




    def get_peaks_newx(self):

        # 's.get_peaks_newx()'
        # 1000 loops, best of 3: 399 usec per loop
        # with reserve: 1000 loops, best of 3: 315 usec per loop


        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef _OpenSwathDataAccessHelper tmp
        cdef shared_ptr[_OSSpectrum] tmp_s

        tmp_s = tmp.convertToSpectrumPtr( deref(spec_) )

        cdef unsigned int n = spec_.size()
        cdef shared_ptr[_OSBinaryDataArray] qq1 = tmp_s.get().getIntensityArray()
        cdef shared_ptr[_OSBinaryDataArray] qq2 = tmp_s.get().getMZArray()

        return qq1.get().data.size()


    def get_peaks_new3(self):
        # 's.get_peaks_new3()'
        # 1000 loops, best of 3: 410 usec per loop

        # 


        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef _OpenSwathDataAccessHelper tmp
        cdef shared_ptr[_OSSpectrum] tmp_s

        tmp_s = tmp.convertToSpectrumPtr( deref(spec_) )

        cdef unsigned int n = spec_.size()
        cdef shared_ptr[_OSBinaryDataArray] qq1 = tmp_s.get().getIntensityArray()
        cdef shared_ptr[_OSBinaryDataArray] qq2 = tmp_s.get().getMZArray()

        ## # We use a memory view to get the data from the raw data
        ## # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        ## # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef double * raw_ptr = address(qq1.get().data[0]) 
        cdef double[:] fda_view = <double[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(fda_view) # numpy array refer to the underlying buffer without copy
        cdef double * raw_ptr2 = address(qq2.get().data[0]) 
        cdef double[:] xx = <double[:n]>raw_ptr2 # cast to memoryview, refer to the underlying buffer without copy
        xarr2 = np.asarray(xx) # numpy array refer to the underlying buffer without copy
        return xarr, xarr2



    def get_peaks_new2(self):


        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()

        cdef libcpp_vector[float] intensities
        cdef libcpp_vector[double] mzs

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        while it != spec_.end():
            mzs.push_back( deref(it).getMZ() )
            intensities.push_back( deref(it).getIntensity() )
            inc(it)

        # return mzs, intensities



        ## # We use a memory view to get the data from the raw data
        ## # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        ## # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef float * raw_ptr = address(intensities[0]) 
        cdef float[:] fda_view = <float[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(fda_view) # numpy array refer to the underlying buffer without copy
        cdef double * raw_ptr2 = address(mzs[0]) 
        cdef double[:] xx = <double[:n]>raw_ptr2 # cast to memoryview, refer to the underlying buffer without copy
        xarr2 = np.asarray(xx) # numpy array refer to the underlying buffer without copy
        return xarr, xarr2

    def get_peaks_new(self):


        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()

        cdef libcpp_vector[float] intensities
        cdef libcpp_vector[double] mzs

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        while it != spec_.end():
            mzs.push_back( deref(it).getMZ() )
            intensities.push_back( deref(it).getIntensity() )
            inc(it)

        return mzs, intensities



        ## cdef _FloatDataArray * fda_ = self.inst.get()
        ## cdef unsigned int n = fda_.size()
        ##  
        ## # Obtain a raw ptr to the beginning of the C++ array
        ## cdef libcpp_vector[float] * vec_ptr = <libcpp_vector[float]*> fda_
        ## cdef float * raw_ptr =  address(deref(vec_ptr)[0]) 

        ## # We use a memory view to get the data from the raw data
        ## # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        ## # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        ## cdef float[:] fda_view = <float[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        ## xarr = np.asarray(fda_view) # numpy array refer to the underlying buffer without copy
        ## return xarr

    def get_peaks(self):

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

    def set_peaks_fast2(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):
        # 1000 loops, best of 3: 1.8 msec per loop

        cdef _OSSpectrum * _s = new _OSSpectrum()
        cdef shared_ptr[_OSSpectrum] tmp_s = shared_ptr[_OSSpectrum](_s)
        cdef shared_ptr[_OSBinaryDataArray] int_arr = tmp_s.get().getIntensityArray()
        cdef shared_ptr[_OSBinaryDataArray] mz_arr = tmp_s.get().getMZArray()

        cdef int N
        N = data_mz.size

        cdef libcpp_vector[double] * int_data = & int_arr.get().data
        cdef libcpp_vector[double] * mz_data = & mz_arr.get().data

        # We use "assign" to directly to copy the numpy array
        cdef double * array_start_mz = <double*>data_mz.data
        mz_data.assign(array_start_mz, array_start_mz + N)
        cdef double * array_start_int = <double*>data_i.data
        int_data.assign(array_start_int, array_start_int + N)

        cdef _MSSpectrum * spec_ = self.inst.get()
        spec_.clear(False)

        cdef _OpenSwathDataAccessHelper tmp
        tmp.convertToOpenMSSpectrum(tmp_s, deref(spec_) )
        spec_.updateRanges()


    def set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):
        # 1000 loops, best of 3: 1.09 msec per loop
        # 1000 loops, best of 3: 781 usec per loop

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


    def set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):
        # 1000 loops, best of 3: 1.09 msec per loop
        # 1000 loops, best of 3: 657 usec per loop

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


    def set_peaks_x(self, peaks):

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"

        mzs, intensities = peaks
        assert len(mzs) == len(intensities), "Input vectors for set_peaks need to have the same length (mz and intensity vector)"

        if isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
                      mzs.dtype == np.float64 and intensities.dtype == np.float32 and \
                      mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :

                        self.set_peaks_fast_df(mzs, intensities)
                        return

        elif isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
                        mzs.dtype == np.float64 and intensities.dtype == np.float64 and \
                        mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :

                        self.set_peaks_fast_dd(mzs, intensities)
                        return

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.clear(0) # empty vector , keep meta data
        spec_.reserve(<int>len(mzs)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float I
        cdef int N
        N = len(mzs)

        for i in range(N):
            mz = mzs[i]
            intensity  = intensities[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()

    def set_peaks(self, peaks):
        # 100 loops, best of 3: 6.22 msec per loop
        # python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand(100000); s=pyopenms.MSSpectrum(); p=pyopenms.Peak1D()' 'for val in d: p.setMZ(val); p.setIntensity(val); s.push_back(p)'
        # 100 loops, best of 3: 17.3 msec per loop

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"

        mzs, intensities = peaks
        assert len(mzs) == len(intensities), "Input vectors for set_peaks need to have the same length (mz and intensity vector)"

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
            intensity  = intensities[i]
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

