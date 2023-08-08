from Types cimport *
from AbsoluteQuantitationMethod cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethodFile:

        AbsoluteQuantitationMethodFile() except + nogil 
        AbsoluteQuantitationMethodFile(AbsoluteQuantitationMethodFile &) except + nogil  # compiler

        void load(const String& filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) except + nogil 
        void store(const String& filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) except + nogil 
