from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from RichPeak2D cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/BaseFeature.h>" namespace "OpenMS":

    cdef cppclass BaseFeature(UniqueIdInterface):
        # wrap-inherits:
        #    UniqueIdInterface

        BaseFeature()  nogil except +
        BaseFeature(BaseFeature &) # wrap-ignore

        Real getQuality()  nogil except +
        void setQuality(Real q) nogil except +

        Real getWidth() nogil except +
        void setWidth(Real q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        bool operator==(BaseFeature) nogil except +
        bool operator!=(BaseFeature) nogil except +

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
        void clearMetaInfo() nogil except +

