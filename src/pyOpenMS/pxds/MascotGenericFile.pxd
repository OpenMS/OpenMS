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
        void store(const String & filename, MSExperiment & experiment) nogil except +
        # NAMESPACE # void store(std::ostream & os, const String & filename, MSExperiment & experiment) nogil except +
        void load(const String & filename, MSExperiment & exp) nogil except +
        libcpp_pair[ String, String ] getHTTPPeakListEnclosure(const String & filename) nogil except +
        void updateMembers_() nogil except +

