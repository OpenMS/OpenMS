from Types cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile "OpenMS::SiriusMSFile":
        SiriusMSFile() nogil except +
        SiriusMSFile(SiriusMSFile) nogil except + #wrap-ignore
        # NAMESPACE # void store(MSExperiment[Peak1D, ChromatogramPeak] & spectra, OpenMS::String & msfile) nogil except +

