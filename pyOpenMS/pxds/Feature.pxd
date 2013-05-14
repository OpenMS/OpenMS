from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DataValue cimport *
from ConvexHull2D cimport *
from PeptideIdentification cimport *

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
        Real getOverallQuality()  nogil except +
        void setOverallQuality(Real q) nogil except +

        Real getWidth() nogil except +
        void setWidth(Real q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +

        libcpp_vector[Feature] getSubordinates() nogil except +
        void setSubordinates(libcpp_vector[Feature]) nogil except +

        ConvexHull2D getConvexHull()                 nogil except +
        libcpp_vector[ConvexHull2D] getConvexHulls() nogil except +
        void setConvexHulls(libcpp_vector[ConvexHull2D]) nogil except +

        # returns a mutable reference to the PeptideIdentification vector
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        # sets the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except +

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
