from Types cimport *
from AbsoluteQuantitationMethod cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethodFile:

        AbsoluteQuantitationMethodFile()  nogil except +
        AbsoluteQuantitationMethodFile(AbsoluteQuantitationMethodFile)  nogil except + #wrap-ignore

        void load(const String& filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) nogil except +
        void store(const String& filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) nogil except +
