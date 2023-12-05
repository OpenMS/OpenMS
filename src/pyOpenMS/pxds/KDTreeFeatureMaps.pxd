from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Feature cimport *
from DefaultParamHandler cimport *
from TransformationModelLowess cimport *
from MSExperiment cimport *
from ConsensusMap cimport *
from FeatureMap cimport *
from BaseFeature cimport *
# from KDTree cimport *
# from KDTreeFeatureNode cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>" namespace "OpenMS":

    cdef cppclass KDTreeFeatureMaps(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        KDTreeFeatureMaps() except + nogil  # wrap-doc:Stores a set of features, together with a 2D tree for fast search
        KDTreeFeatureMaps(libcpp_vector[ FeatureMap ] & maps, Param & param) except + nogil 
        KDTreeFeatureMaps(libcpp_vector[ ConsensusMap ] & maps, Param & param) except + nogil 
        void addMaps(libcpp_vector[ FeatureMap ] & maps) except + nogil  # wrap-doc:Add `maps` and balance kd-tree
        void addMaps(libcpp_vector[ ConsensusMap ] & maps) except + nogil 
        # POINTER # void addFeature(Size mt_map_index, BaseFeature * feature) except + nogil 
        # POINTER # BaseFeature * feature(Size i) except + nogil 
        double rt(Size i) except + nogil 
        double mz(Size i) except + nogil 
        float intensity(Size i) except + nogil 
        Int charge(Size i) except + nogil 
        Size mapIndex(Size i) except + nogil 
        Size size() except + nogil 
        Size treeSize() except + nogil 
        Size numMaps() except + nogil 
        void clear() except + nogil 
        void optimizeTree() except + nogil 
        void getNeighborhood(Size index,
                             libcpp_vector[ size_t ] & result_indices,
                             double rt_tol,
                             double mz_tol,
                             bool mz_ppm,
                             bool include_features_from_same_map,
                             double max_pairwise_log_fc) except + nogil  # wrap-doc:Fill `result` with indices of all features compatible (wrt. RT, m/z, map index) to the feature with `index`
        void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, libcpp_vector[ size_t ] & result_indices, Size ignored_map_index) except + nogil 
        # void applyTransformations(libcpp_vector[ TransformationModelLowess * ] & trafos) except + nogil 

