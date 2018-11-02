cimport numpy as np
import numpy as np




    def get_peaks(self):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        cdef unsigned int n = chrom_.size()
        cdef np.ndarray[np.float64_t, ndim=1] rts
        rts = np.zeros( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.zeros( (n,), dtype=np.float32)
        cdef _ChromatogramPeak p

        cdef libcpp_vector[_ChromatogramPeak].iterator it = chrom_.begin()
        cdef int i = 0
        while it != chrom_.end():
            rts[i] = deref(it).getRT()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return rts, intensities

    def set_peaks(self, peaks):

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"

        rts, intensities = peaks
        assert len(rts) == len(intensities), "Input vectors for set_peaks need to have the same length (rt and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float32 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(rts, intensities)
        elif isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float64 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(rts, intensities)
        else:
            self._set_peaks_orig(rts, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.clear(0) # empty vector , keep meta data
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef double intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.clear(0) # empty vector , keep meta data
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_orig(self, rts, intensities):


        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.clear(0) # empty vector , keep meta data
        chrom_.reserve(<int>len(rts)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(rts)

        for i in range(N):
            rt = rts[i]
            intensity = intensities[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()

