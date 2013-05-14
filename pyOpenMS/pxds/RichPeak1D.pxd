
from Types cimport *
from Peak1D cimport *
from MetaInfoInterface cimport *

cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

    cdef cppclass RichPeak1D(Peak1D, MetaInfoInterface):
        # wrap-inherits:
        #    Peak1D
        #    MetaInfoInterface

        RichPeak1D() nogil except +
        RichPeak1D(RichPeak1D) nogil except +
        bool operator==(RichPeak1D) nogil except +
        bool operator!=(RichPeak1D) nogil except +

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
