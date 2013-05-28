from Types cimport *
from MSExperiment cimport *
from ProgressLogger cimport *
# from DRange cimport *
from PeakWidthEstimator cimport *
from SILACFilter cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>" namespace "OpenMS":
    
    cdef cppclass SILACFiltering(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        SILACFiltering(SILACFiltering) nogil except + #wrap-ignore
        libcpp_vector[SILACFilter] filters_
        # Result peak_width # PeakWidthEstimator # const
        # STL Construct # std::multimap[ DoubleReal, BlacklistEntry ] blacklist
        SILACFiltering(MSExperiment[ Peak1D, ChromatogramPeak] & exp, Result & , DoubleReal intensity_cutoff, String debug_filebase_) nogil except +
        void addFilter(SILACFilter & filter_) nogil except +
        void filterDataPoints() nogil except +

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>" namespace "OpenMS::SILACFiltering":
    
    cdef cppclass BlacklistEntry "OpenMS::SILACFiltering::BlacklistEntry":
        BlacklistEntry(BlacklistEntry) nogil except + #wrap-ignore
        # DRange[ 2 ] range
        Int charge
        libcpp_vector[ double ] mass_separations
        DoubleReal relative_peak_position

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>" namespace "OpenMS::SILACFiltering":
    
    cdef cppclass SpectrumInterpolation "OpenMS::SILACFiltering::SpectrumInterpolation":
        SpectrumInterpolation() nogil except + # wrap-ignore
        SpectrumInterpolation(SpectrumInterpolation) nogil except + #wrap-ignore
        SpectrumInterpolation(MSSpectrum[Peak1D] & , SILACFiltering & ) nogil except +
        # DoubleReal operator()(DoubleReal mz) nogil except +

