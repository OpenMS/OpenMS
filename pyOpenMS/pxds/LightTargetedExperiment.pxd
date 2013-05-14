from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>" namespace "OpenSwath":

    cdef cppclass LightTransition:
        LightTransition()
        LightTransition(LightTransition)
        int getProductChargeState()
        int charge
        libcpp_string transition_name
        libcpp_string peptide_ref
        double library_intensity
        double product_mz

    cdef cppclass LightModification:
        LightModification()
        LightModification(LightModification)
        int location
        libcpp_string unimod_id

    cdef cppclass LightPeptide:
        LightPeptide()
        LightPeptide(LightPeptide)
        double rt
        int charge
        libcpp_string sequence
        libcpp_string protein_ref
        libcpp_string id
        libcpp_vector[LightModification] modifications

    cdef cppclass LightProtein:
        LightProtein()
        LightProtein(LightProtein)
        libcpp_string id
        libcpp_string sequence

    cdef cppclass LightTargetedExperiment:

        LightTargetedExperiment()
        LightTargetedExperiment(LightTargetedExperiment &)

        libcpp_vector[LightTransition] transitions
        libcpp_vector[LightPeptide] peptides
        libcpp_vector[LightProtein] proteins
        libcpp_vector[LightTransition] getTransitions() 


