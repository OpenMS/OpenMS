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

        void pickExperiment(MSExperiment & chromatograms,
                            FeatureMap & output,
                            TargetedExperiment & transition_exp_,
                            TransformationDescription trafo,
                            MSExperiment & swath_map) nogil except +

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

        void prepareProteinPeptideMaps_(LightTargetedExperiment& transition_exp) nogil except +

        # # void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, LightTargetedExperiment transition_exp,
        # #                                    TransitionGroupMapType& transition_group_map,
        # #                                    TransformationDescription trafo, double rt_extraction_window);

