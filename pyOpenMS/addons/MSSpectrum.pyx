# the following empty line is important

    def get_peaks(self):

        cdef _MSSpectrum[_Peak1D] * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float32_t, ndim=2] peaks
        peaks = np.zeros( [n,2], dtype=np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            peaks[i,0] = deref(it).getMZ()
            peaks[i,1] = deref(it).getIntensity()
            inc(it)
            i += 1

        return peaks

    def set_peaks(self, np.ndarray[np.float32_t, ndim=2] peaks):

        cdef _MSSpectrum[_Peak1D] * spec_ = self.inst.get()

        #cdef int delete_meta = 0
        spec_.clear(0) # delete_meta) # emtpy vector , keep meta data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float I
        cdef int N
        N = peaks.shape[0]


        for i in range(N):
            mz = peaks[i,0]
            I  = peaks[i,1]
            p.setMZ(mz)
            p.setIntensity(<float>I)
            spec_.push_back(p)

        spec_.updateRanges()

    def intensityInRange(self, float mzmin, float mzmax):

        cdef int n
        cdef double I

        cdef _MSSpectrum[_Peak1D] * spec_ = self.inst.get()
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

import numpy as np
