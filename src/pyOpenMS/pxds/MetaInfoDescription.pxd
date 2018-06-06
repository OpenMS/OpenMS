from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from String cimport *
from Types cimport *
from MetaInfoRegistry cimport *
from DataProcessing cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/METADATA/MetaInfoDescription.h>" namespace "OpenMS":

    cdef cppclass MetaInfoDescription(MetaInfoInterface):
        # wrap-inherits:
        #  MetaInfoInterface

        MetaInfoDescription() nogil except +
        MetaInfoDescription(MetaInfoDescription) nogil except +

        bool operator==(MetaInfoDescription) nogil except +
        bool operator!=(MetaInfoDescription) nogil except +

        String getName() nogil except +
        void setName(String name) nogil except +

        libcpp_vector[ shared_ptr[DataProcessing] ] getDataProcessing() nogil except +
        void setDataProcessing(libcpp_vector[ shared_ptr[DataProcessing] ]) nogil except +

        # cython has a problem with inheritance of overloaded methods,
        # so we do not declare them here, but separately in each derived
        # class which we want to be wrapped:
        #
        #void getKeys(libcpp_vector[String] & keys)
        #void getKeys(libcpp_vector[unsigned int] & keys)
        #DataValue getMetaValue(unsigned int) nogil except +
        #DataValue getMetaValue(String) nogil except +
        #void setMetaValue(unsigned int, DataValue) nogil except +
        #void setMetaValue(String, DataValue) nogil except +
        #bool metaValueExists(String) nogil except +
        #bool metaValueExists(unsigned int) nogil except +
        #void removeMetaValue(String) nogil except +
        #void removeMetaValue(unsigned int) nogil except +
