cimport numpy as np
import numpy as np


    def get_peaks(self):

        cdef _MSChromatogram[_ChromatogramPeak] * chrom_ = self.inst.get()

        cdef unsigned int n = chrom_.size()
        cdef np.ndarray[np.float32_t, ndim=2] peaks
        peaks = np.zeros( [n,2], dtype=np.float32)
        cdef _ChromatogramPeak p

        cdef libcpp_vector[_ChromatogramPeak].iterator it = chrom_.begin()
        cdef int i = 0
        while it != chrom_.end():
            peaks[i,0] = deref(it).getRT()
            peaks[i,1] = deref(it).getIntensity()
            inc(it)
            i += 1

        return peaks

    def set_peaks(self, np.ndarray[np.float32_t, ndim=2] peaks):

        cdef _MSChromatogram[_ChromatogramPeak] * chrom_ = self.inst.get()

        chrom_.clear(0) # emtpy vector , keep meta data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double mz
        cdef float I
        cdef int N
        N = peaks.shape[0]


        for i in range(N):
            mz = peaks[i,0]
            I  = peaks[i,1]
            p.setRT(mz)
            p.setIntensity(<float>I)
            chrom_.push_back(p)

        chrom_.updateRanges()

