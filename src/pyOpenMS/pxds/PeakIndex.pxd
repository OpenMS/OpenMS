from Types cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/KERNEL/PeakIndex.h>" namespace "OpenMS":
    
    cdef cppclass PeakIndex "OpenMS::PeakIndex":
        # wrap-doc:
                #  Index of a peak or feature
                #  
                #  This struct can be used to store both peak or feature indices
                
        PeakIndex() except + nogil 
        PeakIndex(PeakIndex &) except + nogil  # compiler
        Size peak
        Size spectrum
        PeakIndex(Size peak) except + nogil 
        PeakIndex(Size spectrum, Size peak) except + nogil 
        bool isValid() except + nogil  # wrap-doc:Returns if the current peak ref is valid
        void clear() except + nogil  # wrap-doc:Invalidates the current index
        Feature getFeature(FeatureMap & map_) except + nogil 
            # wrap-doc:
                #  Returns the feature (or consensus feature) corresponding to this index
                #  
                #  This method is intended for arrays of features e.g. FeatureMap
                #  
                #  The main advantage of using this method instead accessing the data directly is that range
                #  check performed in debug mode
                #  
                #  :raises:
                #    Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)

        Peak1D getPeak(MSExperiment & map_) except + nogil 
            # wrap-doc:
                #  Returns a peak corresponding to this index
                #  
                #  This method is intended for arrays of DSpectra e.g. MSExperiment
                #  
                #  The main advantage of using this method instead accessing the data directly is that range
                #  check performed in debug mode
                #  
                #  :raises:
                #    Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)

        MSSpectrum getSpectrum(MSExperiment & map_) except + nogil 
            # wrap-doc:
                #  Returns a spectrum corresponding to this index
                #  
                #  This method is intended for arrays of DSpectra e.g. MSExperiment
                #  
                #  The main advantage of using this method instead accessing the data directly is that range
                #  check performed in debug mode
                #  
                #  :raises:
                #    Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)

        bool operator==(PeakIndex & rhs) except + nogil 
        bool operator!=(PeakIndex & rhs) except + nogil 

