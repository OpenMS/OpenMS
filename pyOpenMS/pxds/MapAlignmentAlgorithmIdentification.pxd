from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from Param cimport *
from Feature cimport *
from FeatureMap cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmIdentification:

        MapAlignmentAlgorithmIdentification() nogil except +
        void alignFeatureMaps(libcpp_vector[FeatureMap[Feature]] & features, libcpp_vector[TransformationDescription] & trafos) nogil except +
        void fitModel(String model_type, Param & params, libcpp_vector[TransformationDescription] & trafos) nogil except +

