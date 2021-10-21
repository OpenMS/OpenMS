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

        OpenSwathHelper() nogil except + # compiler
        OpenSwathHelper(OpenSwathHelper &) nogil except + # compiler

        bool checkSwathMapAndSelectTransitions(
                                MSExperiment & exp, 
                                TargetedExperiment & targeted_exp,
                                TargetedExperiment & transition_exp_used,
                                double min_upper_edge_dist
                                ) nogil except +


        libcpp_pair[double, double] estimateRTRange(LightTargetedExperiment exp) nogil except +
            # wrap-doc:
                #   Computes the min and max retention time value
                #   -----
                #   Estimate the retention time span of a targeted experiment by returning the min/max values in retention time as a pair
                #   -----
                #   :returns: A std `pair` that contains (min,max)

        String computePrecursorId(const String & transition_group_id, int isotope) nogil except +
            # wrap-doc:
                #   Computes unique precursor identifier
                #   -----
                #   Uses transition_group_id and isotope number to compute a unique precursor
                #   id of the form "groupID_Precursor_ix" where x is the isotope number, e.g.
                #   the monoisotopic precursor would become "groupID_Precursor_i0"
                #   -----
                #   :param transition_group_id: Unique id of the transition group (peptide/compound)
                #   :param isotope: Precursor isotope number
                #   :returns: Unique precursor identifier

        # static std::map<std::string, double> simpleFindBestFeature(
        #     OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
        #     bool useQualCutoff = false, double qualCutoff = 0.0);


