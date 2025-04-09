from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from Feature cimport *
from OpenSwathScoring cimport *

cdef extern from "<OpenMS/KERNEL/MRMFeature.h>" namespace "OpenMS":

    cdef cppclass MRMFeature(Feature):
        #
        # wrap-inherits:
        #   Feature

        MRMFeature() except + nogil 
        MRMFeature(MRMFeature &) except + nogil 

        OpenSwath_Scores getScores() except + nogil  # wrap-doc:Returns all peakgroup scores
        void setScores(OpenSwath_Scores s) except + nogil  # wrap-doc:Sets all peakgroup scores

        Feature getFeature(String key) except + nogil  # wrap-doc:Returns a specified feature
        void addFeature(Feature & f, String key) except + nogil  # wrap-doc:Adds an feature from a single chromatogram into the feature
        libcpp_vector[Feature] getFeatures() except + nogil  # wrap-doc:Returns all the features
        void getFeatureIDs(libcpp_vector[String] & result) except + nogil  # wrap-doc:Returns a list of IDs of available features

        Feature getPrecursorFeature(String key) except + nogil  # wrap-doc:Returns a specified precursor feature
        void addPrecursorFeature(Feature & f, String key) except + nogil  # wrap-doc:Adds a precursor feature from a single chromatogram into the feature
        void getPrecursorFeatureIDs(libcpp_vector[String] & result) except + nogil  # wrap-doc:Returns a list of IDs of available precursor features

        bool operator==(MRMFeature) except + nogil 
        bool operator!=(MRMFeature) except + nogil 
