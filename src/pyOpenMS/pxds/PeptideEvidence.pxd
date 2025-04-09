from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/PeptideEvidence.h>" namespace "OpenMS":
    cdef cppclass PeptideEvidence :

        PeptideEvidence() except + nogil 
        PeptideEvidence(PeptideEvidence &) except + nogil 

        # const members
        ## int UNKNOWN_POSITION
        ## int N_TERMINAL_POSITION
        ## char UNKNOWN_AA
        ## char N_TERMINAL_AA
        ## char C_TERMINAL_AA

        void setStart(Int start) except + nogil  # wrap-doc:Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with `N_TERMINAL_AA`
        Int getStart() except + nogil  # wrap-doc:Returns the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        void setEnd(Int end) except + nogil  # wrap-doc:Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA
        Int getEnd() except + nogil  # wrap-doc:Returns the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        void setAABefore(char rhs) except + nogil  # wrap-doc:Sets the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA
        char getAABefore() except + nogil  # wrap-doc:Returns the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned
        void setAAAfter(char rhs) except + nogil  # wrap-doc:Sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA
        char getAAAfter() except + nogil  # wrap-doc:Returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned
        void setProteinAccession(String s) except + nogil  # wrap-doc:Sets the protein accession the peptide matches to. If not available set to empty string
        String getProteinAccession() except + nogil  # wrap-doc:Returns the protein accession the peptide matches to. If not available the empty string is returned

        bool hasValidLimits() except + nogil  # wrap-doc:Start and end numbers in evidence represent actual numeric indices

        bool operator==(PeptideEvidence & rhs) except + nogil 
        bool operator!=(PeptideEvidence & rhs) except + nogil 

