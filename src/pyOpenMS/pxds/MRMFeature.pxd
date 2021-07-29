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

        OpenSwath_Scores getScores() nogil except + Get all peakgroup scores
        void setScores(OpenSwath_Scores s) nogil except + Set all peakgroup scores

        Feature getFeature(String key) nogil except + Get a specified feature
        void addFeature(Feature & f, String key) nogil except + Adds an feature from a single chromatogram into the feature
        libcpp_vector[Feature] getFeatures() nogil except +
        void getFeatureIDs(libcpp_vector[String] & result) nogil except + Get a list of IDs of available features

        Feature getPrecursorFeature(String key) nogil except + Get a specified precursor feature
        void addPrecursorFeature(Feature & f, String key) nogil except + Adds a precursor feature from a single chromatogram into the feature
        void getPrecursorFeatureIDs(libcpp_vector[String] & result) nogil except + Get a list of IDs of available precursor features

        bool operator==(MRMFeature) nogil except +
        bool operator!=(MRMFeature) nogil except +

