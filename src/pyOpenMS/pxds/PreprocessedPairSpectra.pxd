from Types cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":

    cdef cppclass OPXL_PreprocessedPairSpectra "OpenMS::OPXLDataStructs::PreprocessedPairSpectra":

        OPXL_PreprocessedPairSpectra(Size size) nogil except +
        OPXL_PreprocessedPairSpectra(OPXL_PreprocessedPairSpectra &) nogil except +

        MSExperiment spectra_linear_peaks
        MSExperiment spectra_xlink_peaks
        MSExperiment spectra_all_peaks

