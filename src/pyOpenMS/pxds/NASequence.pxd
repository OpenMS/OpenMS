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
        NASequence(NASequence) nogil except + # wrap-ignore

        bool operator==(const NASequence & rhs) nogil except +
        bool operator!=(const NASequence & rhs) nogil except +
        bool operator<(const NASequence & rhs) nogil except +

        libcpp_vector[ const Ribonucleotide * ] getSequence() nogil except +

        #Ribonucleotide * operator[](size_t index) nogil except +



        # check if sequence is empty
        bool empty() nogil except +

        void setSequence(const libcpp_vector[ const Ribonucleotide *] & seq) nogil except +

        # returns the peptide as string with modifications embedded in brackets
        String toString() nogil except +

        # sets the 5' modification
        void setFivePrimeMod(const RibonucleotideChainEnd * modification) nogil except +

        # returns the name (ID) of the N-terminal modification, or an empty string if none is set
        const RibonucleotideChainEnd * getFivePrimeMod() nogil except +

        # sets the 3' modification
        void setThreePrimeMod(const RibonucleotideChainEnd * modification) nogil except +

        const RibonucleotideChainEnd * getThreePrimeMod() nogil except +

        # returns the residue at position index
        const Ribonucleotide * get(Size index) nogil except +

        # set the residue at position index
        void set(size_t index, const Ribonucleotide * r) nogil except +

        # returns the formula of the peptide
        EmpiricalFormula getFormula() nogil except + # wrap-doc:convenience function with ResidueType=Full and charge = 0 by default
        EmpiricalFormula getFormula(NASFragmentType type_, Int charge) nogil except +

        # returns the average weight of the peptide
        double getAverageWeight() nogil except +
        double getAverageWeight(NASFragmentType type_, Int charge) nogil except +

        # returns the mono isotopic weight of the peptide
        double getMonoWeight() nogil except +
        double getMonoWeight(NASFragmentType type_, Int charge) nogil except +

        # returns the number of residues
        Size size() nogil except +

        # returns a peptide sequence of the first index residues
        NASequence getPrefix(Size length) nogil except +

        # returns a peptide sequence of the last index residues
        NASequence getSuffix(Size length) nogil except +

        # returns a peptide sequence of number residues, beginning at position index
        NASequence getSubsequence(Size start, Size length) nogil except +


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

