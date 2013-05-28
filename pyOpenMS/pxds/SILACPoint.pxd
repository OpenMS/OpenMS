from Types cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACPoint.h>" namespace "OpenMS":
    
    cdef cppclass SILACPoint "OpenMS::SILACPoint":
        SILACPoint() nogil except +
        SILACPoint(SILACPoint) nogil except + #wrap-ignore
        DoubleReal mz
        DoubleReal rt
        # libcpp_vector[ libcpp_vector[ double ] ] mz_positions
        # libcpp_vector[ libcpp_vector[ double ] ] intensities
        libcpp_vector[ double ] mass_shifts
        Int charge
        Int isotopes_per_peptide
        DoubleReal quality

