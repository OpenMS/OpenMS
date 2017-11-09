from Types cimport *
from PeptideHit cimport *
from MSExperiment cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    
    cdef cppclass XLPrecursor "OpenMS::OPXLDataStructs::XLPrecursor":
        XLPrecursor(XLPrecursor) nogil except + #wrap-ignore
        float precursor_mass
        unsigned int alpha_index
        unsigned int beta_index

