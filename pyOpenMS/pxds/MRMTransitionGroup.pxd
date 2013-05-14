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

        Size size() nogil except +
        String getTransitionGroupID() nogil except +
        void setTransitionGroupID(String id) nogil except +

        libcpp_vector[TransitionT] getTransitions() nogil except +
        void addTransition(TransitionT & transition, String key) nogil except +
        TransitionT getTransition(String key) nogil except +

