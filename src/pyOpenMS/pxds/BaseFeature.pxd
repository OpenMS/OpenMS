from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from RichPeak2D cimport *
from UniqueIdInterface cimport *
from BaseFeature cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/KERNEL/BaseFeature.h>" namespace "OpenMS":

    cdef cppclass BaseFeature(UniqueIdInterface):
        # wrap-inherits:
        #    UniqueIdInterface

        BaseFeature()  nogil except +
        BaseFeature(BaseFeature &) nogil except +

        float getQuality()  nogil except +
        void setQuality(float q) nogil except +

        float getWidth() nogil except +
        void setWidth(float q) nogil except +

        Int getCharge() nogil except +
        void setCharge(Int q) nogil except +
        AnnotationState getAnnotationState() nogil except +

        # returns a mutable reference to the PeptideIdentification vector
        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except +
        # sets the PeptideIdentification vector
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except +

        bool operator==(BaseFeature) nogil except +
        bool operator!=(BaseFeature) nogil except +

        void getKeys(libcpp_vector[String] & keys) nogil except +
        void getKeys(libcpp_vector[unsigned int] & keys) nogil except + # wrap-as:getKeysAsIntegers
        DataValue getMetaValue(unsigned int) nogil except +
        DataValue getMetaValue(String) nogil except +
        void setMetaValue(unsigned int, DataValue) nogil except +
        void setMetaValue(String, DataValue) nogil except +
        bool metaValueExists(String) nogil except +
        bool metaValueExists(unsigned int) nogil except +
        void removeMetaValue(String) nogil except +
        void removeMetaValue(unsigned int) nogil except +
        void clearMetaInfo() nogil except +

cdef extern from "<OpenMS/KERNEL/BaseFeature.h>" namespace "OpenMS::BaseFeature":
    
    cdef enum AnnotationState "OpenMS::BaseFeature::AnnotationState":
        FEATURE_ID_NONE
        FEATURE_ID_SINGLE
        FEATURE_ID_MULTIPLE_SAME
        FEATURE_ID_MULTIPLE_DIVERGENT
        SIZE_OF_ANNOTATIONSTATE
