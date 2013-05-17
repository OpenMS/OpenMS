from Types cimport *
from libcpp cimport bool
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentEvaluationAlgorithm "OpenMS::MapAlignmentEvaluationAlgorithm":
        # wrap-ignore
        # ABSTRACT class
        MapAlignmentEvaluationAlgorithm() nogil except +
        MapAlignmentEvaluationAlgorithm(MapAlignmentEvaluationAlgorithm) nogil except + #wrap-ignore
        # NAMESPACE # void evaluate(ConsensusMap & conensus_map_in, ConsensusMap & consensus_map_gt, DoubleReal & rt_dev, DoubleReal & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge, DoubleReal & out) nogil except +
        # NAMESPACE # bool isSameHandle(FeatureHandle & lhs, FeatureHandle & rhs, DoubleReal & rt_dev, DoubleReal & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge) nogil except +
        void registerChildren() nogil except +

