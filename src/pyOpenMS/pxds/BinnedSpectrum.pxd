from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *
# from SparseVector cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>" namespace "OpenMS":

    cdef cppclass BinnedSpectrum:

        BinnedSpectrum() nogil except +
        BinnedSpectrum(BinnedSpectrum) nogil except +
        BinnedSpectrum(MSSpectrum, float size, bool unit_ppm, UInt spread) nogil except +

        bool operator==(BinnedSpectrum & rhs) nogil except +
        bool operator!=(BinnedSpectrum & rhs) nogil except +

        float getBinSize() nogil except +
        UInt getBinSpread() nogil except +
        UInt getBinIndex(float mz) nogil except +
        float getBinLowerMZ(size_t i) nogil except +

        vector[Precursor] getPrecursors() nogil except +

        bool isCompatible(BinnedSpectrum & a, BinnedSpectrum & b) nogil except +
