from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from RichPeak2D cimport *
from UniqueIdInterface cimport *
from BaseFeature cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/KERNEL/BaseFeature.h>" namespace "OpenMS":

    cdef cppclass BaseFeature(UniqueIdInterface, RichPeak2D):
        # wrap-inherits:
        #    UniqueIdInterface
        #    RichPeak2D

        BaseFeature()  nogil except +
        BaseFeature(BaseFeature &) nogil except +

        float getQuality()  nogil except + # wrap-doc:Returns the overall quality
        void setQuality(float q) nogil except + # wrap-doc:Sets the overall quality

        float getWidth() nogil except + # wrap-doc:Returns the features width (full width at half max, FWHM)
        void setWidth(float q) nogil except + # wrap-doc:Sets the width of the feature (FWHM)

        Int getCharge() nogil except + # wrap-doc:Returns the charge state
        void setCharge(Int q) nogil except + # wrap-doc:Sets the charge state
        AnnotationState getAnnotationState() nogil except + # wrap-doc:State of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only

        libcpp_vector[PeptideIdentification] getPeptideIdentifications() nogil except + # wrap-doc:Returns the PeptideIdentification vector
        
        void setPeptideIdentifications(libcpp_vector[PeptideIdentification] & peptides) nogil except + # wrap-doc:Sets the PeptideIdentification vector

        bool operator==(BaseFeature) nogil except +
        bool operator!=(BaseFeature) nogil except +

cdef extern from "<OpenMS/KERNEL/BaseFeature.h>" namespace "OpenMS::BaseFeature":
    
    cdef enum AnnotationState "OpenMS::BaseFeature::AnnotationState":
        FEATURE_ID_NONE
        FEATURE_ID_SINGLE
        FEATURE_ID_MULTIPLE_SAME
        FEATURE_ID_MULTIPLE_DIVERGENT
        SIZE_OF_ANNOTATIONSTATE
