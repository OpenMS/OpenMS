from Types cimport *
from libcpp cimport bool
from MapAlignmentEvaluationAlgorithm cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentEvaluationAlgorithmRecall(MapAlignmentEvaluationAlgorithm) :
        # wrap-inherits:
        #  MapAlignmentEvaluationAlgorithm
        MapAlignmentEvaluationAlgorithmRecall() except + nogil 
        # private
        MapAlignmentEvaluationAlgorithmRecall(MapAlignmentEvaluationAlgorithmRecall) except + nogil  #wrap-ignore

        # NAMESPACE # void evaluate(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_gt, double & rt_dev, double & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge, double & out) except + nogil 
        # POINTER # MapAlignmentEvaluationAlgorithm * create() except + nogil 
       