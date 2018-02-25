from Types cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from String cimport *
from StringList cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/SiriusMSConverter.h>" namespace "OpenMS":
    
    cdef cppclass SiriusMSFile "OpenMS::SiriusMSFile":
        SiriusMSFile() nogil except +
        SiriusMSFile(SiriusMSFile) nogil except + #wrap-ignore
        # void store(MSExperiment spectra, String msfile, libcpp_map[size_t, StringList] map_precursor_to_adducts) nogil except +

