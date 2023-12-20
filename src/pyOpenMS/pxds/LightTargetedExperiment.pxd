from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string

cdef extern from "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>" namespace "OpenSwath":

    cdef cppclass LightTransition:
        LightTransition() except + nogil 
        LightTransition(LightTransition &) except + nogil 
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

        int getProductChargeState() except + nogil 
        bool isProductChargeStateSet() except + nogil 
        libcpp_string getNativeID() except + nogil 
        libcpp_string getPeptideRef() except + nogil 
        double getLibraryIntensity() except + nogil 
        void setLibraryIntensity(double l) except + nogil 
        double getProductMZ() except + nogil 
        double getPrecursorMZ() except + nogil 

        libcpp_string getCompoundRef() except + nogil 

        # Detecting / quantifying / identifying transitions
        void setDetectingTransition (bool d) except + nogil 
        bool isDetectingTransition() except + nogil 
        void setQuantifyingTransition (bool q) except + nogil 
        bool isQuantifyingTransition() except + nogil 
        void setIdentifyingTransition (bool i) except + nogil 
        bool isIdentifyingTransition() except + nogil 

    cdef cppclass LightModification:
        LightModification() except + nogil 
        LightModification(LightModification &) except + nogil 

        int location
        int unimod_id

    cdef cppclass LightCompound:
        LightCompound() except + nogil 
        LightCompound(LightCompound &) except + nogil 

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

        void setDriftTime(double d) except + nogil 
        double getDriftTime() except + nogil 

        int getChargeState() except + nogil 
        bool isPeptide() except + nogil 
        void setChargeState(int ch) except + nogil 

    cdef cppclass LightProtein:
        LightProtein() except + nogil 
        LightProtein(LightProtein &) except + nogil 
        libcpp_string id
        libcpp_string sequence

    cdef cppclass LightTargetedExperiment:

        LightTargetedExperiment() except + nogil 
        LightTargetedExperiment(LightTargetedExperiment &) except + nogil 

        libcpp_vector[LightTransition] transitions
        libcpp_vector[LightCompound] compounds
        libcpp_vector[LightProtein] proteins
        libcpp_vector[LightTransition] getTransitions()  except + nogil 

        libcpp_vector[LightCompound] getCompounds() except + nogil 
        libcpp_vector[LightProtein] getProteins() except + nogil 

        LightCompound getCompoundByRef(libcpp_string & ref) except + nogil 
        LightCompound getPeptideByRef(libcpp_string & ref) except + nogil 

