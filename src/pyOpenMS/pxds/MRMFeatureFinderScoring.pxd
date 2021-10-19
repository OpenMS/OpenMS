from Types cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureMap cimport *
from SwathMap cimport *
from TargetedExperiment cimport *
from LightTargetedExperiment cimport *
from TransformationDescription cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from SpectrumAccessOpenMS cimport *
from ISpectrumAccess cimport *
from MRMTransitionGroup cimport *

from SpectrumAccessOpenMS cimport *
from SpectrumAccessOpenMSCached cimport *
###  TODO: enable once we have these classes
### # from SpectrumAccessOpenMSInMemory cimport *
### # from SpectrumAccessQuadMZTransforming cimport *

# typedef OpenSwath::LightTransition TransitionType;
# typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFinderScoring(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MRMFeatureFinderScoring() nogil except +
        # copy constructor of 'MRMFeatureFinderScoring' is implicitly deleted because field 'diascoring_' has an inaccessible copy constructor
        MRMFeatureFinderScoring(MRMFeatureFinderScoring &) nogil except + # wrap-ignore

        void pickExperiment(MSExperiment & chromatograms,
                            FeatureMap & output,
                            TargetedExperiment & transition_exp_,
                            TransformationDescription trafo,
                            MSExperiment & swath_map) nogil except +
            # wrap-doc:
                #   Pick features in one experiment containing chromatogram
                #   -----
                #   Function for for wrapping in Python, only uses OpenMS datastructures and does not return the map
                #   -----
                #   :param chromatograms: The input chromatograms
                #   :param output: The output features with corresponding scores
                #   :param transition_exp: The transition list describing the experiment
                #   :param trafo: Optional transformation of the experimental retention time to the normalized retention time space used in the transition list
                #   :param swath_map: Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted

        ## void pickExperiment(shared_ptr[ SpectrumAccessOpenMS ] input_chrom,
        ##                     FeatureMap& output,
        ##                     LightTargetedExperiment& transition_exp,
        ##                     TransformationDescription trafo,
        ##                     SwathMap swath_map,
        ##                     libcpp_map[String, MRMTransitionGroup[MSChromatogram, LightTransition] ] & transition_group_map)

        void setStrictFlag(bool flag) nogil except +

        # void setMS1Map( shared_ptr[ ISpectrumAccess ] ms1_map) nogil except +
        void setMS1Map( shared_ptr[ SpectrumAccessOpenMS ] ms1_map) nogil except +
        void setMS1Map( shared_ptr[ SpectrumAccessOpenMSCached ] ms1_map) nogil except +
        ### # void setMS1Map( shared_ptr[ SpectrumAccessOpenMSInMemory ] ms1_map) nogil except +
        ### # void setMS1Map( shared_ptr[ SpectrumAccessQuadMZTransforming ] ms1_map) nogil except +

        void scorePeakgroups(MRMTransitionGroup[MSChromatogram, LightTransition] transition_group,
                             TransformationDescription trafo,
                             libcpp_vector[ SwathMap ] swath_maps,
                             FeatureMap& output,
                             bool ms1only) nogil except +
            # wrap-doc:
                #   Score all peak groups of a transition group
                #   -----
                #   Iterate through all features found along the chromatograms of the transition group and score each one individually
                #   -----
                #   :param transition_group: The MRMTransitionGroup to be scored (input)
                #   :param trafo: Optional transformation of the experimental retention time
                #               to the normalized retention time space used in the
                #               transition list
                #   :param swath_maps: Optional SWATH-MS (DIA) map corresponding from which
                #                    the chromatograms were extracted. Use empty map if no
                #                    data is available
                #   :param output: The output features with corresponding scores (the found
                #                features will be added to this FeatureMap)
                #   :param ms1only: Whether to only do MS1 scoring and skip all MS2 scoring

        void prepareProteinPeptideMaps_(LightTargetedExperiment& transition_exp) nogil except +
            # wrap-doc:
                #   Prepares the internal mappings of peptides and proteins
                #   -----
                #   Calling this method _is_ required before calling scorePeakgroups
                #   -----
                #   :param transition_exp: The transition list describing the experiment

        # # void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, LightTargetedExperiment transition_exp,
        # #                                    TransitionGroupMapType& transition_group_map,
        # #                                    TransformationDescription trafo, double rt_extraction_window);

