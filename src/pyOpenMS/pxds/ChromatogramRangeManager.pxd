from Types cimport *
from RangeManager cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramRangeManager.h>" namespace "OpenMS":
    
    cdef cppclass ChromatogramRangeManager:
        # wrap-doc:
        #  Range manager for chromatograms
        #  
        #  This class manages retention time, m/z, and intensity ranges for multiple chromatograms.
        #  It extends the basic RangeManager to provide specialized functionality for chromatogram data.
        #  
        #  The template parameters for the base RangeManager are ordered differently than in SpectrumRangeManager:
        #  - RangeRT (retention time) is the first parameter, as it's the primary dimension for chromatograms
        #  - RangeIntensity is the second parameter
        #  - RangeMZ is the third parameter
        
        ChromatogramRangeManager() except + nogil
        ChromatogramRangeManager(ChromatogramRangeManager &) except + nogil
        
        void clearRanges() except + nogil
        
        # Range accessors
        double getMinRT() except + nogil
        double getMaxRT() except + nogil
        double getMinMZ() except + nogil
        double getMaxMZ() except + nogil
        double getMinIntensity() except + nogil
        double getMaxIntensity() except + nogil