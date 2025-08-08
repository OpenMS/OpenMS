from Types cimport *
from RangeManager cimport *
from MSSpectrum cimport *
from libcpp.set cimport set as libcpp_set

cdef extern from "<OpenMS/KERNEL/SpectrumRangeManager.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumRangeManager:
        # wrap-doc:
        #  Advanced range manager for MS spectra with separate ranges for each MS level
        #  
        #  This class extends the basic RangeManager to provide separate range tracking for different MS levels
        #  (MS1, MS2, etc.). It manages four types of ranges:
        #  - m/z (mass-to-charge ratio)
        #  - intensity
        #  - retention time (RT)
        #  - ion mobility
        #  
        #  A global range is tracked for all MS levels, and additional ranges are maintained for each specific MS level.
        #  This allows for efficient querying of ranges for specific MS levels, which is useful for visualization,
        #  filtering, and processing operations that need to work with specific MS levels.
        
        SpectrumRangeManager() except + nogil
        SpectrumRangeManager(SpectrumRangeManager &) except + nogil
        
        bool operator==(SpectrumRangeManager &) except + nogil
        bool operator!=(SpectrumRangeManager &) except + nogil
        
        void clearRanges() except + nogil
        libcpp_set[UInt] getMSLevels() except + nogil
        void extendRT(double rt, UInt ms_level) except + nogil
        void extendMZ(double mz, UInt ms_level) except + nogil
        void extendUnsafe(const MSSpectrum& spectrum, UInt ms_level) except + nogil
        
        # Range accessors
        double getMinRT() except + nogil
        double getMaxRT() except + nogil
        double getMinMZ() except + nogil
        double getMaxMZ() except + nogil
        double getMinIntensity() except + nogil
        double getMaxIntensity() except + nogil
        double getMinMobility() except + nogil
        double getMaxMobility() except + nogil