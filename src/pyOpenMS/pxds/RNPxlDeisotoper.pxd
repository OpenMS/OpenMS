from libcpp cimport bool
from Types cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlDeisotoper.h>" namespace "OpenMS":

    cdef cppclass RNPxlDeisotoper:
        RNPxlDeisotoper() nogil except +
        RNPxlDeisotoper(RNPxlDeisotoper) nogil except + # wrap-ignore

# COMMENT: wrap static methods
cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlDeisotoper.h>" namespace "OpenMS::RNPxlDeisotoper":
        # static members
        void deisotopeAndSingleCharge(MSSpectrum& in, double fragment_tolerance, bool fragment_unit_ppm, int min_charge, int max_charge, bool keep_only_deisotoped, unsigned int min_isopeaks,  unsigned int max_isopeaks, bool make_single_charged = true, bool annotate_charge = false) nogil except + # wrap-attach:RNPxlDeisotoper

