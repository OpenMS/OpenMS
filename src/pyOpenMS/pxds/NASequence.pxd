from libcpp cimport bool
from Types cimport *
from String cimport *
from Residue cimport *
from EmpiricalFormula cimport *
from Ribonucleotide cimport *
from Map cimport *

cdef extern from "<OpenMS/CHEMISTRY/NASequence.h>" namespace "OpenMS":

    ctypedef Ribonucleotide RibonucleotideChainEnd

    cdef cppclass NASequence "OpenMS::NASequence":
        # wrap-hash:
        #   toString().c_str()
        #
        # wrap-doc:
        #   Representation of an RNA sequence
        #   This class represents nucleic acid sequences in OpenMS. An NASequence
        #   instance primarily contains a sequence of ribonucleotides. 

        NASequence() nogil except +
        NASequence(NASequence &) nogil except +

        bool operator==(const NASequence & rhs) nogil except +
        bool operator!=(const NASequence & rhs) nogil except +
        bool operator<(const NASequence & rhs) nogil except +

        libcpp_vector[ const Ribonucleotide * ] getSequence() nogil except +

        const Ribonucleotide * operator[](size_t index) nogil except +

        bool empty() nogil except + # wrap-doc:Check if sequence is empty

        void setSequence(const libcpp_vector[ const Ribonucleotide *] & seq) nogil except +

        String toString() nogil except + # wrap-doc:Returns the peptide as string with modifications embedded in brackets

        void setFivePrimeMod(const RibonucleotideChainEnd * modification) nogil except + # wrap-doc:Sets the 5' modification

        const RibonucleotideChainEnd * getFivePrimeMod() nogil except + # wrap-doc:Returns the name (ID) of the N-terminal modification, or an empty string if none is set

        void setThreePrimeMod(const RibonucleotideChainEnd * modification) nogil except + # wrap-doc:Sets the 3' modification

        const RibonucleotideChainEnd * getThreePrimeMod() nogil except +

        const Ribonucleotide * get(Size index) nogil except + # wrap-doc:Returns the residue at position index

        void set(size_t index, const Ribonucleotide * r) nogil except + # wrap-doc:Sets the residue at position index

        EmpiricalFormula getFormula() nogil except + # wrap-doc:Returns the formula of the peptide
        EmpiricalFormula getFormula(NASFragmentType type_, Int charge) nogil except +

        double getAverageWeight() nogil except + # wrap-doc:Returns the average weight of the peptide
        double getAverageWeight(NASFragmentType type_, Int charge) nogil except +

        double getMonoWeight() nogil except + # wrap-doc:Returns the mono isotopic weight of the peptide
        double getMonoWeight(NASFragmentType type_, Int charge) nogil except +

        Size size() nogil except + # wrap-doc:Returns the number of residues

        NASequence getPrefix(Size length) nogil except + # wrap-doc:Returns a peptide sequence of the first index residues

        NASequence getSuffix(Size length) nogil except + # wrap-doc:Returns a peptide sequence of the last index residues

        NASequence getSubsequence(Size start, Size length) nogil except + # wrap-doc:Returns a peptide sequence of number residues, beginning at position index


# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/NASequence.h>" namespace "OpenMS::NASequence":
        
        
        # static members
        NASequence fromString(const String & s) nogil except +  # wrap-attach:NASequence

cdef extern from "<OpenMS/CHEMISTRY/NASequence.h>" namespace "OpenMS::NASequence":
    cdef enum NASFragmentType "OpenMS::NASequence::NASFragmentType":
        #wrap-attach:
        #    NASequence
        Full
        Internal
        FivePrime
        ThreePrime
        AIon
        BIon
        CIon
        XIon
        YIon
        ZIon
        Precursor
        BIonMinusH20
        YIonMinusH20
        BIonMinusNH3
        YIonMinusNH3
        NonIdentified
        Unannotated
        WIon
        AminusB
        DIon
        SizeOfNASFragmentType

