from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *

cdef extern from "<OpenMS/KERNEL/BinnedSpectrum.h>" namespace "OpenMS":

    cdef cppclass BinnedSpectrum:
        BinnedSpectrum() except + nogil 
        BinnedSpectrum(BinnedSpectrum &) except + nogil 
        BinnedSpectrum(MSSpectrum, float size, bool unit_ppm, UInt spread, float offset) except + nogil 

        bool operator==(BinnedSpectrum & rhs) except + nogil 
        bool operator!=(BinnedSpectrum & rhs) except + nogil 

        float getBinSize() except + nogil  # wrap-doc:Returns the bin size
        UInt getBinSpread() except + nogil  # wrap-doc:Returns the bin spread
        UInt getBinIndex(float mz) except + nogil  # wrap-doc:Returns the bin index of a given m/z position
        float getBinLowerMZ(size_t i) except + nogil  # wrap-doc:Returns the lower m/z of a bin given its index
        float getBinIntensity(double mz) except + nogil  # wrap-doc:Returns the bin intensity at a given m/z position

        libcpp_vector[Precursor] getPrecursors() except + nogil 

        bool isCompatible(BinnedSpectrum & a, BinnedSpectrum & b) except + nogil 
       
        float getOffset() except + nogil  # wrap-doc:Returns offset
