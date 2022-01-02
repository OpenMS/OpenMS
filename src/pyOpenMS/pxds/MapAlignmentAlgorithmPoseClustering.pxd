from ChromatogramPeak cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmPoseClustering(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MapAlignmentAlgorithmPoseClustering() nogil except +

        # private
        MapAlignmentAlgorithmPoseClustering(MapAlignmentAlgorithmPoseClustering &) nogil except + # wrap-ignore

        void align(FeatureMap,
                   TransformationDescription &
                   ) nogil except +

        void align(MSExperiment,
                   TransformationDescription &
                   ) nogil except +

        void setReference (FeatureMap) nogil except + # wrap-doc:Sets the reference for the alignment
        void setReference (MSExperiment) nogil except +


