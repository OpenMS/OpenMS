from Types cimport *
from OPXLDataStructs cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":

    cdef cppclass XLPrecursor "OpenMS::OPXLDataStructs::XLPrecursor":

        XLPrecursor() except + nogil 
        XLPrecursor(XLPrecursor &) except + nogil 


        float precursor_mass
        unsigned int alpha_index
        unsigned int beta_index
