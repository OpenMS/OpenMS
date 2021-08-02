from Types cimport *
from libcpp cimport bool
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>" namespace "OpenMS":
    
    cdef cppclass MapAlignmentEvaluationAlgorithm "OpenMS::MapAlignmentEvaluationAlgorithm":
        # wrap-ignore
        # ABSTRACT class
        # no-pxd-import
        MapAlignmentEvaluationAlgorithm() nogil except +
        # private
        MapAlignmentEvaluationAlgorithm(MapAlignmentEvaluationAlgorithm) nogil except + #wrap-ignore

        # NAMESPACE # void evaluate(ConsensusMap & conensus_map_in, ConsensusMap & consensus_map_gt, double & rt_dev, double & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge, double & out) nogil except +
        # NAMESPACE # bool isSameHandle(FeatureHandle & lhs, FeatureHandle & rhs, double & rt_dev, double & mz_dev, Peak2D::IntensityType & int_dev, bool use_charge) nogil except +
        void registerChildren() nogil except + # wrap-doc:Register all derived classes in this method
