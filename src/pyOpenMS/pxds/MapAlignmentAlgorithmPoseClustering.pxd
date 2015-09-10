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

        void align(FeatureMap,
                   TransformationDescription &
                   ) nogil except +

        void align(MSExperiment[Peak1D,ChromatogramPeak],
                   TransformationDescription &
                   ) nogil except +

        void setReference (FeatureMap) nogil except +
        void setReference (MSExperiment[Peak1D,ChromatogramPeak]) nogil except +


