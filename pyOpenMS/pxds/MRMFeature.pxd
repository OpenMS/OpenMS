from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from Feature cimport *

cdef extern from "<OpenMS/KERNEL/MRMFeature.h>" namespace "OpenMS":

    cdef cppclass MRMFeature:

        MRMFeature() nogil except +
        MRMFeature(MRMFeature) nogil except +

        # libcpp_map[String, double] getScores() nogil except +
        double getScore(String name) nogil except +
        void addScore(String name, double score) nogil except +

        Feature getFeature(String key) nogil except +
        void addFeature(Feature & f, String key) nogil except +
        libcpp_vector[Feature] getFeatures() nogil except +
        void getFeatureIDs(libcpp_vector[String] & result) nogil except +

