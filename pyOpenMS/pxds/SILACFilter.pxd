from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from SILACFiltering cimport *
from SILACPattern cimport *
# from IsotopeDistributionCache cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>" namespace "OpenMS":
    
    cdef cppclass SILACFilter "OpenMS::SILACFilter":
        SILACFilter(SILACFilter) nogil except + #wrap-ignore
        SILACFilter(libcpp_vector[ double ] mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide, DoubleReal intensity_cutoff, DoubleReal intensity_correlation, bool allow_missing_peaks) nogil except +
        libcpp_vector[ double ] getPeakPositions() nogil except +
        libcpp_vector[ double ]  getExpectedMzShifts() nogil except + 
        libcpp_vector[ SILACPattern ]  getElements() nogil except +
        Int getCharge() nogil except +
        libcpp_vector[ double ]  getMassSeparations() nogil except +

