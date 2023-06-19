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
        KDTreeFeatureMaps() nogil except + # wrap-doc:Stores a set of features, together with a 2D tree for fast search
        KDTreeFeatureMaps(libcpp_vector[ FeatureMap ] & maps, Param & param) nogil except +
        KDTreeFeatureMaps(libcpp_vector[ ConsensusMap ] & maps, Param & param) nogil except +
        void addMaps(libcpp_vector[ FeatureMap ] & maps) nogil except + # wrap-doc:Add `maps` and balance kd-tree
        void addMaps(libcpp_vector[ ConsensusMap ] & maps) nogil except +
        # POINTER # void addFeature(Size mt_map_index, BaseFeature * feature) nogil except +
        # POINTER # BaseFeature * feature(Size i) nogil except +
        double rt(Size i) nogil except +
        double mz(Size i) nogil except +
        float intensity(Size i) nogil except +
        Int charge(Size i) nogil except +
        Size mapIndex(Size i) nogil except +
        Size size() nogil except +
        Size treeSize() nogil except +
        Size numMaps() nogil except +
        void clear() nogil except +
        void optimizeTree() nogil except +
        void getNeighborhood(Size index,
                             libcpp_vector[ size_t ] & result_indices,
                             double rt_tol,
                             double mz_tol,
                             bool mz_ppm,
                             bool include_features_from_same_map,
                             double max_pairwise_log_fc) nogil except + # wrap-doc:Fill `result` with indices of all features compatible (wrt. RT, m/z, map index) to the feature with `index`
        void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, libcpp_vector[ size_t ] & result_indices, Size ignored_map_index) nogil except +
        # void applyTransformations(libcpp_vector[ TransformationModelLowess * ] & trafos) nogil except +

