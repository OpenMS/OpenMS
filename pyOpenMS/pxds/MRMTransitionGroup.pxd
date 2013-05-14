from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from MRMFeature cimport *
from ReactionMonitoringTransition cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/KERNEL/MRMTransitionGroup.h>" namespace "OpenMS":

    cdef cppclass MRMTransitionGroup[SpectrumT, TransitionT]:

        # wrap-instances:
        #   MRMTransitionGroup := MRMTransitionGroup[MSSpectrum[Peak1D], ReactionMonitoringTransition]

        MRMTransitionGroup() nogil except +
        MRMTransitionGroup(MRMTransitionGroup) nogil except +

        Size size() nogil except+

        String getTransitionGroupID() nogil except+

        void setTransitionGroupID(String tr_gr_id)

        libcpp_vector[TransitionT] getTransitions() nogil except+

        libcpp_vector[TransitionT] getTransitionsMuteable()

        void addTransition(TransitionT transition, String key)

        TransitionT getTransition(String key) 

        bool hasTransition(String key)

        libcpp_vector[SpectrumT] getChromatograms() nogil except+

        void addChromatogram(SpectrumT chromatogram, String key)

        SpectrumT getChromatogram(String key)

        bool hasChromatogram(String key)

        libcpp_vector[MRMFeature] getFeatures() nogil except+

        libcpp_vector[MRMFeature] getFeaturesMuteable()

        void addFeature(MRMFeature feature)

        void getLibraryIntensity(libcpp_vector[double] result) nogil except+

"""
        Size size() nogil except +
        String getTransitionGroupID() nogil except +
        void setTransitionGroupID(String id) nogil except +

        libcpp_vector[TransitionT] getTransitions() nogil except +
        void addTransition(TransitionT transition, String key) nogil except +
        TransitionT getTransition(String key) nogil except +
        """

