from Types cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/SqMassFile.h>" namespace "OpenMS":
    
    cdef cppclass SqMassFile "OpenMS::SqMassFile":
        SqMassFile() nogil except +
        SqMassFile(SqMassFile) nogil except + #wrap-ignore
        void load(String & filename, MSExperiment & map_) nogil except +
        void store(String & filename, MSExperiment & map_) nogil except +

