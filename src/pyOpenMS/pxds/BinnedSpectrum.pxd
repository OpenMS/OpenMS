from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>" namespace "OpenMS":

    cdef cppclass BinnedSpectrum:
        BinnedSpectrum() nogil except + #wrap-ignore
        BinnedSpectrum(BinnedSpectrum &) nogil except + #wrap-ignore
        BinnedSpectrum(MSSpectrum, float size, bool unit_ppm, UInt spread, float offset) nogil except +

        bool operator==(BinnedSpectrum & rhs) nogil except +
        bool operator!=(BinnedSpectrum & rhs) nogil except +

        float getBinSize() nogil except +
        UInt getBinSpread() nogil except +
        UInt getBinIndex(float mz) nogil except +
        float getBinLowerMZ(size_t i) nogil except +
        float getBinIntensity(double mz) nogil except +

        libcpp_vector[Precursor] getPrecursors() nogil except +

        bool isCompatible(BinnedSpectrum & a, BinnedSpectrum & b) nogil except +
