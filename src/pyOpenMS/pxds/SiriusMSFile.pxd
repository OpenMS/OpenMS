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
        # void store(PeakMap spectra, String msfile, libcpp_map[size_t, BaseFeature] feature_ms2_spectra_map, bool feature_only, int isotope_pattern_iterations, bool no_mt_info) nogil except +

