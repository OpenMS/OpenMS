from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FeatureMap cimport *
from TransformationDescription cimport *
from BinaryTreeNode cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmTreeGuided(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger
        MapAlignmentAlgorithmTreeGuided() except + nogil
        MapAlignmentAlgorithmTreeGuided(MapAlignmentAlgorithmTreeGuided&) except + nogil
        void treeGuidedAlignment(const libcpp_vector[BinaryTreeNode]& tree,
                                 libcpp_vector[FeatureMap]& feature_maps_transformed,
                                 libcpp_vector[libcpp_vector[double]]& maps_ranges,
                                 FeatureMap& map_transformed,
                                 libcpp_vector[Size]& trafo_order) except + nogil
        void align(libcpp_vector[FeatureMap]& data,
                   libcpp_vector[TransformationDescription]& transformations) except + nogil
        void computeTrafosByOriginalRT(libcpp_vector[FeatureMap]& feature_maps,
                                       FeatureMap& map_transformed,
                                       libcpp_vector[TransformationDescription]& transformations,
                                       const libcpp_vector[Size]& trafo_order) except + nogil
        @staticmethod
        void buildTree(libcpp_vector[FeatureMap]& feature_maps,
                       libcpp_vector[BinaryTreeNode]& tree,
                       libcpp_vector[libcpp_vector[double]]& maps_ranges) except + nogil
        @staticmethod
        void computeTransformedFeatureMaps(libcpp_vector[FeatureMap]& feature_maps,
                                           const libcpp_vector[TransformationDescription]& transformations) except + nogil
