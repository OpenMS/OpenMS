from Types cimport *
from libcpp cimport bool
from MapAlignmentEvaluationAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentEvaluationAlgorithmRecall(MapAlignmentEvaluationAlgorithm) :
        # wrap-inherits:
        #  MapAlignmentEvaluationAlgorithm
        MapAlignmentEvaluationAlgorithmRecall() nogil except +
        MapAlignmentEvaluationAlgorithmRecall(MapAlignmentEvaluationAlgorithmRecall) nogil except + #wrap-ignore
        # NAMESPACE # void evaluate(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_gt, DoubleReal & rt_dev, DoubleReal & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge, DoubleReal & out) nogil except +
        # POINTER # MapAlignmentEvaluationAlgorithm * create() nogil except +
        String getProductName() nogil except +

