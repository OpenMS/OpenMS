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
        MetaInfoDescription(MetaInfoDescription& ) nogil except +

        bool operator==(MetaInfoDescription) nogil except +
        bool operator!=(MetaInfoDescription) nogil except +

        String getName() nogil except + # wrap-doc:Returns the name of the peak annotations
        void setName(String name) nogil except + # wrap-doc:Sets the name of the peak annotations

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except + # wrap-doc:Returns a reference to the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) nogil except + # wrap-doc:Sets the description of the applied processing

