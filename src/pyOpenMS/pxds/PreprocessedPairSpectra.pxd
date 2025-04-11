from Types cimport *
from MSRun cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":

    cdef cppclass OPXL_PreprocessedPairSpectra "OpenMS::OPXLDataStructs::PreprocessedPairSpectra":

        OPXL_PreprocessedPairSpectra(Size size) except + nogil 
        OPXL_PreprocessedPairSpectra(OPXL_PreprocessedPairSpectra &) except + nogil 

        MSRun spectra_linear_peaks
        MSRun spectra_xlink_peaks
        MSRun spectra_all_peaks

