from Types cimport *
from MSExperiment cimport *
from IMSDataConsumer cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS":
    
    cdef cppclass SqMassFile "OpenMS::SqMassFile":
        SqMassFile() nogil except +
        SqMassFile(SqMassFile) nogil except + #wrap-ignore
        void load(const String & filename, MSExperiment & map_) nogil except +
        void store(const String & filename, MSExperiment & map_) nogil except +
        # TODO !
        void transform(const String&,
                       IMSDataConsumer[Peak1D, ChromatogramPeak] *,
                       bool skip_full_count,
                       bool skip_first_pass) nogil except + # wrap-ignore

