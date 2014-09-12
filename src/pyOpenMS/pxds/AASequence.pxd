from libcpp cimport bool
from Types cimport *
from String cimport *
from Residue cimport *
from EmpiricalFormula cimport *
from Map cimport *

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS":

    cdef cppclass AASequence:

        AASequence() nogil except +
        AASequence(AASequence) nogil except + # wrap-ignore


        AASequence operator+(AASequence)    nogil except +
        AASequence iadd(AASequence)   nogil except + # wrap-as:operator+=

        # check if sequence is empty
        bool empty() nogil except +

        # returns the peptide as string with modifications embedded in brackets
        String toString() nogil except +

        # returns the peptide as string without any modifications
        String toUnmodifiedString() nogil except +

        # set the modification of the residue at position index
        void setModification(Size index, String modification) nogil except +

        # sets the N-terminal modification
        void setNTerminalModification(String modification) nogil except +

        # returns the Id of the N-term modification nogil except + an empty string is returned if none was set
        String getNTerminalModification() nogil except +

        # sets the C-terminal modification
        void setCTerminalModification(String modification) nogil except +

        # returns the Id of the C-term modification nogil except + an empty string is returned if none was set
        String getCTerminalModification() nogil except +

        # returns a pointer to the residue, which is at position index
        Residue getResidue(SignedSize index) nogil except +

        # returns a pointer to the residue, which is at position index
        Residue getResidue(Size index) nogil except +

        # returns the formula of the peptide
        EmpiricalFormula getFormula(ResidueType type_, Int charge) nogil except +

        # returns the average weight of the peptide
        double getAverageWeight(ResidueType type_, Int charge) nogil except +

        # returns the mono isotopic weight of the peptide
        double getMonoWeight(ResidueType type_, Int charge) nogil except +

        # returns the number of residues
        Size size() nogil except +

        # returns a peptide sequence of the first index residues
        AASequence getPrefix(Size index) nogil except +

        # returns a peptide sequence of the last index residues
        AASequence getSuffix(Size index) nogil except +

        # returns a peptide sequence of number residues, beginning at position index
        AASequence getSubsequence(Size index, UInt number) nogil except +

        # compute frequency table of amino acids
        void getAAFrequencies(Map[String, size_t]) nogil except + # wrap-ignore

        # returns true if the peptude contains the given residue
        bool has(Residue residue) nogil except +

        # returns true if the peptide contains the given peptide
        # @note c-term and n-term mods are ignored
        bool hasSubsequence(AASequence peptide) nogil except +

        # returns true if the peptide has the given prefix
        # n-term mod is also checked (c-term as well, if prefix is of same length)
        bool hasPrefix(AASequence peptide) nogil except +

        # returns true if the peptide has the given suffix
        # c-term mod is also checked (n-term as well, if suffix is of same length)
        bool hasSuffix(AASequence peptide) nogil except +

        # predicate which is true if the peptide is N-term modified
        bool hasNTerminalModification() nogil except +

        # predicate which is true if the peptide is C-term modified
        bool hasCTerminalModification() nogil except +

        # returns true if any of the residues is modified
        bool isModified() nogil except +

        # returns true if the residue at the position is modified
        bool isModified(Size index) nogil except +

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS::AASequence":

        # static members
        AASequence fromString(String s, bool permissive) nogil except +  # wrap-attach:AASequence

