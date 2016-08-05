from Types cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureMap cimport *
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
# typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; 

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>" namespace "OpenMS":

    cdef cppclass MRMFeatureFinderScoring(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ProgressLogger

        MRMFeatureFinderScoring() nogil except +

        void pickExperiment(MSExperiment[Peak1D, ChromatogramPeak] & chromatograms,
                            FeatureMap & output,
                            TargetedExperiment & transition_exp_,
                            TransformationDescription trafo,
                            MSExperiment[Peak1D, ChromatogramPeak] & swath_map) nogil except +

        # void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap& output, OpenSwath::LightTargetedExperiment& transition_exp,
        #                     TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType& transition_group_map);

        # void pickExperiment(shared_ptr[ ISpectrumAccess ] ms1_map inp,
        #                     FeatureMap& output,
        #                     LightTargetedExperiment& transition_exp,
        #                     TransformationDescription trafo, 
        #                     OpenSwath::SpectrumAccessPtr swath_map, 
        #                     TransitionGroupMapType& transition_group_map)

        void setStrictFlag(bool flag) nogil except +

        # void setMS1Map( shared_ptr[ ISpectrumAccess ] ms1_map) nogil except +
        void setMS1Map( shared_ptr[ SpectrumAccessOpenMS ] ms1_map) nogil except +
        void setMS1Map( shared_ptr[ SpectrumAccessOpenMSCached ] ms1_map) nogil except +
        ### # void setMS1Map( shared_ptr[ SpectrumAccessOpenMSInMemory ] ms1_map) nogil except +
        ### # void setMS1Map( shared_ptr[ SpectrumAccessQuadMZTransforming ] ms1_map) nogil except +

        # void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo,
        #                      shared_ptr[ ISpectrumAccess ] swath_map, FeatureMap& output, bool ms1only) nogil except +
        void scorePeakgroups(MRMTransitionGroup[MSSpectrum[ChromatogramPeak], LightTransition] transition_group,
                             TransformationDescription trafo,
                             shared_ptr[ SpectrumAccessOpenMS ] swath_map,
                             FeatureMap& output, bool ms1only) nogil except +
        void scorePeakgroups(MRMTransitionGroup[MSSpectrum[ChromatogramPeak], LightTransition] transition_group,
                             TransformationDescription trafo,
                             shared_ptr[ SpectrumAccessOpenMSCached ] swath_map,
                             FeatureMap& output, bool ms1only) nogil except +
        ### # void scorePeakgroups(MRMTransitionGroup[MSSpectrum[ChromatogramPeak], LightTransition] transition_group,
        ### #                      TransformationDescription trafo,
        ### #                      shared_ptr[ SpectrumAccessOpenMSInMemory ] swath_map,
        ### #                      FeatureMap& output, bool ms1only) nogil except +
        ### # void scorePeakgroups(MRMTransitionGroup[MSSpectrum[ChromatogramPeak], LightTransition] transition_group,
        ### #                      TransformationDescription trafo,
        ### #                      shared_ptr[ SpectrumAccessQuadMZTransforming ] swath_map,
        ### #                      FeatureMap& output, bool ms1only) nogil except +
          
        void prepareProteinPeptideMaps_(LightTargetedExperiment& transition_exp) nogil except +

