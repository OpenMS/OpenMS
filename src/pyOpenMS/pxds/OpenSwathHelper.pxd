from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *
from TransformationDescription cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>" namespace "OpenMS":

    cdef cppclass OpenSwathHelper:

        bool checkSwathMapAndSelectTransitions(
                                MSExperiment & exp, 
                                TargetedExperiment & targeted_exp,
                                TargetedExperiment & transition_exp_used,
                                double min_upper_edge_dist
                                ) nogil except +


        libcpp_pair[double, double] estimateRTRange(LightTargetedExperiment exp) nogil except +

        String computePrecursorId(const String & transition_group_id, int isotope) nogil except +

        # static std::map<std::string, double> simpleFindBestFeature(
        #     OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
        #     bool useQualCutoff = false, double qualCutoff = 0.0);


