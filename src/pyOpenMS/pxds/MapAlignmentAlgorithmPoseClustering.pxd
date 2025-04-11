from ChromatogramPeak cimport *
from DefaultParamHandler cimport *
from Feature cimport *
from FeatureMap cimport *
from MSRun cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TransformationDescription cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>" namespace "OpenMS":

    cdef cppclass MapAlignmentAlgorithmPoseClustering(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger

        MapAlignmentAlgorithmPoseClustering() except + nogil 

        # private
        MapAlignmentAlgorithmPoseClustering(MapAlignmentAlgorithmPoseClustering &) except + nogil  # wrap-ignore

        void align(FeatureMap,
                   TransformationDescription &
                   ) except + nogil 

        void align(MSRun,
                   TransformationDescription &
                   ) except + nogil 

        void setReference (FeatureMap) except + nogil  # wrap-doc:Sets the reference for the alignment
        void setReference (MSRun) except + nogil 


