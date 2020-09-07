from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoInterface cimport *
from DataProcessing cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoDescription.h>" namespace "OpenMS":

    cdef cppclass MetaInfoDescription(MetaInfoInterface):
        # wrap-inherits:
        #   MetaInfoInterface

        MetaInfoDescription() nogil except +
        MetaInfoDescription(MetaInfoDescription) nogil except +

        bool operator==(MetaInfoDescription) nogil except +
        bool operator!=(MetaInfoDescription) nogil except +

        String getName() nogil except +
        void setName(String name) nogil except +

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) nogil except +

