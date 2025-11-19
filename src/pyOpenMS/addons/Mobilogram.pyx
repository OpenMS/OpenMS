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

