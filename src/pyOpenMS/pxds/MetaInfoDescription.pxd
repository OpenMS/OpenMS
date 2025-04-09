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
        #  MetaInfoInterface

        MetaInfoDescription() except + nogil 
        MetaInfoDescription(MetaInfoDescription& ) except + nogil 

        bool operator==(MetaInfoDescription) except + nogil 
        bool operator!=(MetaInfoDescription) except + nogil 

        String getName() except + nogil  # wrap-doc:Returns the name of the peak annotations
        void setName(String name) except + nogil  # wrap-doc:Sets the name of the peak annotations

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() except + nogil  # wrap-doc:Returns a reference to the description of the applied processing
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) except + nogil  # wrap-doc:Sets the description of the applied processing

