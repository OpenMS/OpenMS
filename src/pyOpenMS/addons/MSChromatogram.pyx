cimport numpy as np
import numpy as np


    def get_peaks(self):

        cdef _MSChromatogram[_ChromatogramPeak] * chrom_ = self.inst.get()

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

        assert isinstance(peaks, (tuple, list))
        assert len(peaks) == 2

        rts, intensities = peaks
        assert len(rts) == len(intensities)

        cdef _MSChromatogram[_ChromatogramPeak] * chrom_ = self.inst.get()

        chrom_.clear(0) # emtpy vector, keep meta data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float I
        cdef int N
        N = len(rts)

        for i in range(N):
            rt = rts[i]
            intensity  = intensities[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


