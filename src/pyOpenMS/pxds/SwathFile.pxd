from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from MSExperiment cimport *
from MzMLFile cimport *
from MzXMLFile cimport *
from SwathMap cimport *
from StringList cimport *
from ExperimentalSettings cimport *

cdef extern from "<OpenMS/FORMAT/SwathFile.h>" namespace "OpenMS":
    
    cdef cppclass SwathFile(ProgressLogger) :
        # wrap-inherits:
        #  ProgressLogger
        SwathFile(SwathFile) nogil except + #wrap-ignore

        libcpp_vector[ SwathMap ] loadSplit(StringList file_list,
                                            String tmp,
                                            shared_ptr[ ExperimentalSettings ] exp_meta,
                                            String readoptions) nogil except +

        libcpp_vector[ SwathMap ] loadMzML(String file_,
                                           String tmp,
                                           shared_ptr[ ExperimentalSettings ] exp_meta,
                                           String readoptions) nogil except +

        libcpp_vector[ SwathMap ] loadMzXML(String file_,
                                            String tmp,
                                            shared_ptr[ ExperimentalSettings ] exp_meta, 
                                            String readoptions) nogil except +

