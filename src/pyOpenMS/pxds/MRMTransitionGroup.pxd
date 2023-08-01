from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_utf8_string # triggers input conversion provider
from MRMFeature cimport *
from ReactionMonitoringTransition cimport *
from LightTargetedExperiment cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/KERNEL/MRMTransitionGroup.h>" namespace "OpenMS":

    cdef cppclass MRMTransitionGroup[SpectrumT, TransitionT]:

        # wrap-instances:
        #  MRMTransitionGroupCP := MRMTransitionGroup[MSChromatogram, ReactionMonitoringTransition]
        #  LightMRMTransitionGroupCP := MRMTransitionGroup[MSChromatogram, LightTransition]

        MRMTransitionGroup() except + nogil # TODO(whole file)
        MRMTransitionGroup(MRMTransitionGroup[SpectrumT, TransitionT] &) except + nogil 

        Size size() nogil except+

        String getTransitionGroupID() nogil except+
        void setTransitionGroupID(String tr_gr_id) except + nogil 

        libcpp_vector[TransitionT] getTransitions() nogil except+
        libcpp_vector[TransitionT] getTransitionsMuteable() except + nogil 
        void addTransition(TransitionT transition, String key) except + nogil 
        TransitionT getTransition(String key)  except + nogil 
        bool hasTransition(String key) except + nogil 

        libcpp_vector[SpectrumT] getChromatograms() nogil except+
        void addChromatogram(SpectrumT chromatogram, String key) except + nogil 
        SpectrumT getChromatogram(String key) except + nogil 
        bool hasChromatogram(String key) except + nogil 

        libcpp_vector[SpectrumT] getPrecursorChromatograms() nogil except+
        void addPrecursorChromatogram(SpectrumT chromatogram, String key) except + nogil 
        SpectrumT getPrecursorChromatogram(String key) except + nogil 
        bool hasPrecursorChromatogram(String key) except + nogil 

        libcpp_vector[MRMFeature] getFeatures() nogil except+
        libcpp_vector[MRMFeature] getFeaturesMuteable() except + nogil 
        void addFeature(MRMFeature feature) except + nogil 
        MRMFeature getBestFeature() except + nogil 

        void getLibraryIntensity(libcpp_vector[double] result) nogil except+
        MRMTransitionGroup[SpectrumT, TransitionT] subset(libcpp_vector[ libcpp_utf8_string ] tr_ids) except + nogil 

        bool isInternallyConsistent() except + nogil 

        bool chromatogramIdsMatch() except + nogil 
