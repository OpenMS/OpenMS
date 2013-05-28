from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from String cimport *
from MSExperiment cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FORMAT/MascotGenericFile.h>" namespace "OpenMS":
    
    cdef cppclass MascotGenericFile(ProgressLogger, DefaultParamHandler) :
        # wrap-inherits:
        #  ProgressLogger
        #  DefaultParamHandler
        MascotGenericFile() nogil except +
        MascotGenericFile(MascotGenericFile) nogil except + #wrap-ignore
        void store(String & filename, MSExperiment[Peak1D, ChromatogramPeak] & experiment) nogil except +
        # NAMESPACE # void store(std::ostream & os, String & filename, MSExperiment[Peak1D, ChromatogramPeak] & experiment) nogil except +
        void load(String & filename, MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +
        libcpp_pair[ String, String ] getHTTPPeakListEnclosure(String & filename) nogil except +

