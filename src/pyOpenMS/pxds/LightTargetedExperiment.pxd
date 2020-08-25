from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>" namespace "OpenSwath":

    cdef cppclass LightTransition:
        LightTransition() nogil except +
        LightTransition(LightTransition) nogil except +
        libcpp_string transition_name
        libcpp_string peptide_ref
        double library_intensity
        double product_mz
        double precursor_mz

        int fragment_charge
        bool decoy 
        bool detecting_transition
        bool quantifying_transition
        bool identifying_transition

        int getProductChargeState() nogil except +
        bool isProductChargeStateSet() nogil except +
        libcpp_string getNativeID() nogil except +
        libcpp_string getPeptideRef() nogil except +
        double getLibraryIntensity() nogil except +
        void setLibraryIntensity(double l) nogil except +
        double getProductMZ() nogil except +
        double getPrecursorMZ() nogil except +

        libcpp_string getCompoundRef() nogil except +

        # Detecting / quantifying / identifying transitions
        void setDetectingTransition (bool d) nogil except +
        bool isDetectingTransition() nogil except +
        void setQuantifyingTransition (bool q) nogil except +
        bool isQuantifyingTransition() nogil except +
        void setIdentifyingTransition (bool i) nogil except +
        bool isIdentifyingTransition() nogil except +

    cdef cppclass LightModification:
        LightModification() nogil except +
        LightModification(LightModification) nogil except +

        int location
        int unimod_id

    cdef cppclass LightCompound:
        LightCompound() nogil except +
        LightCompound(LightCompound) nogil except +

        double rt
        double drift_time
        int charge
        libcpp_string sequence
        libcpp_vector[libcpp_string] protein_refs
        libcpp_string peptide_group_label
        libcpp_string id
        libcpp_string sum_formula
        libcpp_string compound_name

        libcpp_vector[LightModification] modifications

        void setDriftTime(double d) nogil except +
        double getDriftTime() nogil except +

        int getChargeState() nogil except +
        bool isPeptide() nogil except +
        void setChargeState(int ch) nogil except +

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

        libcpp_vector[LightCompound] getCompounds() nogil except +
        libcpp_vector[LightProtein] getProteins() nogil except +

        LightCompound getCompoundByRef(libcpp_string & ref) nogil except +
        LightCompound getPeptideByRef(libcpp_string & ref) nogil except +

