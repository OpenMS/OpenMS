from Types cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from String cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile "OpenMS::SiriusMSFile":
        SiriusMSFile() nogil except +
        SiriusMSFile(SiriusMSFile) nogil except + #wrap-ignore
        void store(MSExperiment & spectra, String & msfile) nogil except +

