from Types cimport *
from Param cimport *
from AbsoluteQuantitationMethod cimport *
cdef extern from "<OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>" namespace "OpenMS":

    cdef cppclass AbsoluteQuantitationMethodFile:
        # wrap-ignore
        # no-pxd-import

        AbsoluteQuantitationMethodFile()  nogil except +

        void load(string filename, libcpp_vector[ dataAbsoluteQuantitationMethod ]& aqm_list) nogil except +

