from Types cimport *
from PeptideHit cimport *
from MSExperiment cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    
    cdef cppclass OPXL_PreprocessedPairSpectra "OpenMS::OPXLDataStructs::PreprocessedPairSpectra":
        OPXL_PreprocessedPairSpectra(OPXL_PreprocessedPairSpectra) nogil except + #wrap-ignore
        MSExperiment spectra_common_peaks
        MSExperiment spectra_xlink_peaks
        MSExperiment spectra_all_peaks
        OPXL_PreprocessedPairSpectra(Size size) nogil except +

