from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>" namespace "OpenMS":

    cdef cppclass BinnedSpectrum:
        BinnedSpectrum() nogil except +
        BinnedSpectrum(BinnedSpectrum &) nogil except +
        BinnedSpectrum(MSSpectrum, float size, bool unit_ppm, UInt spread, float offset) nogil except +

        bool operator==(BinnedSpectrum & rhs) nogil except +
        bool operator!=(BinnedSpectrum & rhs) nogil except +

        float getBinSize() nogil except + # wrap-doc:Returns the bin size
        UInt getBinSpread() nogil except + # wrap-doc:Returns the bin spread
        UInt getBinIndex(float mz) nogil except + # wrap-doc:Returns the bin index of a given m/z position
        float getBinLowerMZ(size_t i) nogil except + # wrap-doc:Returns the lower m/z of a bin given its index
        float getBinIntensity(double mz) nogil except + # wrap-doc:Returns the bin intensity at a given m/z position

        libcpp_vector[Precursor] getPrecursors() nogil except +

        bool isCompatible(BinnedSpectrum & a, BinnedSpectrum & b) nogil except +
       
        float getOffset() nogil except + # wrap-doc:Returns offset
