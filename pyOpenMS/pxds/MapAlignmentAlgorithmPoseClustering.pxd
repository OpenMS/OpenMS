from Peak1D cimport *
from ChromatogramPeak cimport *
from TransformationDescription cimport *
from MSExperiment cimport *
from String cimport *
from Param cimport *
from Feature cimport *
from FeatureMap cimport *
from ProgressLogger cimport *
from DefaultParamHandler cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmPoseClustering(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MapAlignmentAlgorithmPoseClustering() nogil except +

        void align(FeatureMap[Feature],
                   TransformationDescription &
                   ) nogil except +

        void align(MSExperiment[Peak1D,ChromatogramPeak],
                   TransformationDescription &
                   ) nogil except +

        void setReference (FeatureMap[Feature]) nogil except +
        void setReference (MSExperiment[Peak1D,ChromatogramPeak]) nogil except +
        String getProductName()


