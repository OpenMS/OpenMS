from Types cimport *
from SwathMap cimport *
from StringList cimport *
from String cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/SwathFile.h>" namespace "OpenMS":

    cdef cppclass SwathFile:

        SwathFile() nogil except +
        SwathFile(SwathFile) nogil except +

        libcpp_vector[ SwathMap ] loadSplit(StringList file_list, String tmp,
            shared_ptr[ ExperimentalSettings] exp_meta, String readoptions) nogil except +
        libcpp_vector[ SwathMap ] loadMzML(String file_, String tmp,
            shared_ptr[ ExperimentalSettings] exp_meta, String readoptions) nogil except +
        libcpp_vector[ SwathMap ] loadMzXML(String file_, String tmp,
            shared_ptr[ ExperimentalSettings] exp_meta, String readoptions) nogil except +

