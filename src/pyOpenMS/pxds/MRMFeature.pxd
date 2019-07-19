from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from Feature cimport *
from OpenSwathScoring cimport *

cdef extern from "<OpenMS/KERNEL/MRMFeature.h>" namespace "OpenMS":

    cdef cppclass MRMFeature(Feature):
        #
        # wrap-inherits:
        #    Feature

        MRMFeature() nogil except +
        MRMFeature(MRMFeature &) nogil except +

        OpenSwath_Scores getScores() nogil except +
        void setScores(OpenSwath_Scores s) nogil except +

        Feature getFeature(String key) nogil except +
        void addFeature(Feature & f, String key) nogil except +
        libcpp_vector[Feature] getFeatures() nogil except +
        void getFeatureIDs(libcpp_vector[String] & result) nogil except +

        Feature getPrecursorFeature(String key) nogil except +
        void addPrecursorFeature(Feature & f, String key) nogil except +
        void getPrecursorFeatureIDs(libcpp_vector[String] & result) nogil except +

        bool operator==(MRMFeature) nogil except +
        bool operator!=(MRMFeature) nogil except +

