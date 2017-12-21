from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *

cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethodFile:

        AbsoluteQuantitationMethodFile()  nogil except +
        AbsoluteQuantitationMethodFile(AbsoluteQuantitationMethodFile)  nogil except + #wrap-ignore

        void load(String filename, libcpp_vector[ AbsoluteQuantitationMethod ]& aqm_list) nogil except +

