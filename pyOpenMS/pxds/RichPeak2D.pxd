
from libcpp cimport bool
from Types cimport *
from Peak2D cimport *
from MetaInfoInterface cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS":

    cdef cppclass RichPeak2D(Peak2D, UniqueIdInterface, MetaInfoInterface):
        # wrap-inherits:
        #    Peak2D
        #    UniqueIdInterface
        #    MetaInfoInterface

        RichPeak2D() nogil except +

        bool operator==(RichPeak2D)
        bool operator!=(RichPeak2D)

        # declare again: cython complains for overloaded methods in base
        # classes
        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except +
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
