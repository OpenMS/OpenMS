from Types cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from SpectrumAccessOpenMS cimport *
from ISpectrumAccess cimport *
from MRMTransitionGroup cimport *

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

        void setStrictFlag(bool flag) nogil except +

        void setMS1Map( shared_ptr[ SpectrumAccessOpenMS ] ms1_map) nogil except +
        # TODO for ISpectrumAccess
        # void setMS1Map( shared_ptr[ ISpectrumAccess ] ms1_map) nogil except +

        # TODO for ISpectrumAccess
        # void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo,
        #                      OpenSwath::SpectrumAccessPtr swath_map, FeatureMap& output);

        # void scorePeakgroups(MRMTransitionGroup[MSSpectrum[Peak1D], LightTransition] transition_group, TransformationDescription trafo,
        #                      shared_ptr[ SpectrumAccessOpenMS ] swath_map, FeatureMap& output) nogil except +
          
        void prepareProteinPeptideMaps_(LightTargetedExperiment& transition_exp) nogil except +
           
