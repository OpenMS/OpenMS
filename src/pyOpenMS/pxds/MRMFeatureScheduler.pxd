from FeatureMap cimport *
from String cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureScheduler.h>" namespace "OpenMS":

    cdef cppclass MRMFS_SelectorParameters "OpenMS::MRMFeatureScheduler::SelectorParameters":
        MRMFS_SelectorParameters()
        MRMFS_SelectorParameters(MRMFS_SelectorParameters &) # no-wrap

        Int nn_threshold
        bool locality_weight
        bool select_transition_group
        Int segment_window_length
        Int segment_step_length
        String variable_type
        double optimal_threshold
        libcpp_map[String,String] score_weights

    cdef cppclass MRMFeatureScheduler:

        MRMFeatureScheduler() nogil except +
        MRMFeatureScheduler(MRMFeatureScheduler &) nogil except +

        void schedule_MRMFeaturesScore(FeatureMap& features, FeatureMap& selected_features) nogil except +
        void schedule_MRMFeaturesQMIP(FeatureMap& features, FeatureMap& selected_features) nogil except +

        void setSchedulerParameters(libcpp_vector[SelectorParameters]& parameters) nogil except +
        libcpp_vector[SelectorParameters]& getSchedulerParameters(void) nogil except +
