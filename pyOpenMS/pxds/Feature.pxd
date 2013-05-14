from libcpp cimport bool
from DataValue cimport *


from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/Feature.h>" namespace "OpenMS":

    cdef cppclass Feature(UniqueIdInterface):
        #
        # wrap-inherits:
        #    UniqueIdInterface

        Feature() nogil except +
        Feature(Feature &) nogil except +

        void setMZ(double)  nogil except +
        void setRT(double)  nogil except +
        void setIntensity(double) nogil except +

        double getMZ() nogil except +
        double getRT() nogil except +
        double getIntensity() nogil except +

        Real getQuality(Size index)  nogil except +
        void setQuality(Size index, Real q) nogil except +

        Real getWidth() nogil except +
        void setWidth(Real q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        bool operator==(Feature) nogil except +
        bool operator!=(Feature) nogil except +

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
