from OPXLDataStructs cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS":

    cdef cppclass XLPrecursor "OpenMS::OPXLDataStructs::XLPrecursor":
        XLPrecursor(XLPrecursor) nogil except +
        #wrap-attach:
        #    OPXLDataStructs

        float precursor_mass
        unsigned int alpha_index
        unsigned int beta_index
