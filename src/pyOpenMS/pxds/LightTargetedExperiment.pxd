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
        double precursor_im

        int fragment_charge
        # Note: decoy, detecting_transition, quantifying_transition, identifying_transition
        # are now packed in flags struct - use getter/setter methods below

        # Additional fields for roundtrip TSV/PQP I/O
        int fragment_nr
        # Note: fragment_type is now an enum - use getFragmentType()/setFragmentType() methods
        # Note: annotation is now reconstructed via getAnnotation() - no direct field access
        libcpp_vector[libcpp_string] peptidoforms

        int getProductChargeState() except + nogil
        bool isProductChargeStateSet() except + nogil
        libcpp_string getNativeID() except + nogil
        libcpp_string getPeptideRef() except + nogil
        double getLibraryIntensity() except + nogil
        void setLibraryIntensity(double l) except + nogil
        double getProductMZ() except + nogil
        double getPrecursorMZ() except + nogil
        bool isPrecursorImSet() except + nogil
        double getPrecursorIM() except + nogil

        libcpp_string getCompoundRef() except + nogil

        # Decoy getter/setter
        bool getDecoy() except + nogil
        void setDecoy(bool d) except + nogil

        # Fragment type getter/setter (string interface for backward compatibility)
        libcpp_string getFragmentType() except + nogil
        void setFragmentType(libcpp_string s) except + nogil

        # Annotation reconstruction (computed from fragment_type, fragment_nr, fragment_charge)
        libcpp_string getAnnotation() except + nogil

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
        libcpp_string gene_name
        libcpp_string id
        libcpp_string sum_formula
        libcpp_string compound_name

        # Additional fields for roundtrip TSV/PQP I/O
        libcpp_string label_type
        libcpp_string smiles
        libcpp_string adducts

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

        # Additional fields for roundtrip TSV/PQP I/O
        libcpp_string uniprot_id

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

