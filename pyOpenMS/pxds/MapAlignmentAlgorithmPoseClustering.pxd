from Peak1D cimport *
from ChromatogramPeak cimport *
from TransformationDescription cimport *
from MSExperiment cimport *
from String cimport *
from Param cimport *
from Feature cimport *
from FeatureMap cimport *
from DefaultParamHandler cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmPoseClustering(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        MapAlignmentAlgorithmPoseClustering() nogil except +
        void align(FeatureMap[Feature], TransformationDescription &) nogil except +
        void align(MSExperiment[Peak1D,ChromatogramPeak], TransformationDescription &) nogil except +

        void setReference (FeatureMap[Feature]) nogil except +
        void setReference (MSExperiment[Peak1D,ChromatogramPeak]) nogil except +

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentTransformer:
        pass

       
cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>" namespace "OpenMS::MapAlignmentTransformer":

    void transformFeatureMaps(libcpp_vector[FeatureMap[Feature]], libcpp_vector[TransformationDescription] &) nogil except + #wrap-attach:MapAlignmentTransformer
    void transformPeakMaps(libcpp_vector[MSExperiment[Peak1D,ChromatogramPeak]], libcpp_vector[TransformationDescription]) nogil except + #wrap-attach:MapAlignmentTransformer

    void transformSingleFeatureMap(FeatureMap[Feature] &, TransformationDescription &) nogil except + #wrap-attach:MapAlignmentTransformer
    void transformSinglePeakMap(MSExperiment[Peak1D,ChromatogramPeak] &, TransformationDescription &) nogil except + #wrap-attach:MapAlignmentTransformer

