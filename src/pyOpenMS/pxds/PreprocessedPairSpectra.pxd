from Types cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":

    cdef cppclass OPXL_PreprocessedPairSpectra "OpenMS::OPXLDataStructs::PreprocessedPairSpectra":

        OPXL_PreprocessedPairSpectra(Size size) except + nogil 
        OPXL_PreprocessedPairSpectra(OPXL_PreprocessedPairSpectra &) except + nogil 

        MSExperiment spectra_linear_peaks
        MSExperiment spectra_xlink_peaks
        MSExperiment spectra_all_peaks

