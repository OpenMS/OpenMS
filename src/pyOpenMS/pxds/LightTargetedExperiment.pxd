from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>" namespace "OpenSwath":

    cdef cppclass LightTransition:
        LightTransition() nogil except +
        LightTransition(LightTransition) nogil except +
        int getProductChargeState() nogil except +
        int charge
        bool decoy 
        libcpp_string transition_name
        libcpp_string peptide_ref
        double library_intensity
        double product_mz
        double precursor_mz

        libcpp_string getNativeID() nogil except +
        libcpp_string getPeptideRef() nogil except +
        double getLibraryIntensity() nogil except +
        void setLibraryIntensity(double l) nogil except +
        double getProductMZ() nogil except +
        double getPrecursorMZ() nogil except +

    cdef cppclass LightModification:
        LightModification() nogil except +
        LightModification(LightModification) nogil except +
        int location
        libcpp_string unimod_id

    cdef cppclass LightCompound:
        LightCompound() nogil except +
        LightCompound(LightCompound) nogil except +
        double rt
        int charge
        libcpp_string sequence
        libcpp_vector[libcpp_string] protein_refs
        libcpp_string peptide_group_label
        libcpp_string id
        libcpp_string sum_formula
        libcpp_string compound_name

        libcpp_vector[LightModification] modifications

        int getChargeState() nogil except +

    cdef cppclass LightProtein:
        LightProtein() nogil except +
        LightProtein(LightProtein) nogil except +
        libcpp_string id
        libcpp_string sequence

    cdef cppclass LightTargetedExperiment:

        LightTargetedExperiment() nogil except +
        LightTargetedExperiment(LightTargetedExperiment &) nogil except +

        libcpp_vector[LightTransition] transitions
        libcpp_vector[LightCompound] compounds
        libcpp_vector[LightProtein] proteins
        libcpp_vector[LightTransition] getTransitions()  nogil except +

        libcpp_vector[ LightCompound ]  getCompounds() nogil except +
        libcpp_vector[ LightProtein ]  getProteins() nogil except +

        LightCompound getCompoundByRef(libcpp_string & ref) nogil except +

