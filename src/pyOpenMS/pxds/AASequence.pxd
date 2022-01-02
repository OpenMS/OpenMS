from libcpp cimport bool
from Types cimport *
from String cimport *
from Residue cimport *
from EmpiricalFormula cimport *
from ResidueModification cimport *
from Map cimport *

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS":

    cdef cppclass AASequence:
        # wrap-hash:
        #   toString().c_str()
        #
        # wrap-doc:
        #   Representation of a peptide/protein sequence
        #   This class represents amino acid sequences in OpenMS. An AASequence
        #   instance primarily contains a sequence of residues. 

        AASequence() nogil except +
        AASequence(AASequence &) nogil except +

        AASequence operator+(AASequence) nogil except +
        AASequence iadd(AASequence) nogil except + # wrap-as:operator+=

        # Note that this is a const-ref, so we cannot easily set residues
        Residue operator[](int) nogil except + # wrap-upper-limit:size()

        # check if sequence is empty
        bool empty() nogil except + # wrap-doc:Check if sequence is empty

        # returns the peptide as string with modifications embedded in brackets
        String toString() nogil except + # wrap-doc:Returns the peptide as string with modifications embedded in brackets

        # returns the peptide as string without any modifications
        String toUnmodifiedString() nogil except + # wrap-doc:Returns the peptide as string without any modifications

        # returns the peptide as string without any modifications
        String toUniModString() nogil except + # wrap-doc:Returns the peptide as string with UniMod-style modifications embedded in brackets

        String toBracketString() nogil except + # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation. Uses integer mass by default
        String toBracketString(bool integer_mass) nogil except + # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation
        String toBracketString(bool integer_mass, bool mass_delta) nogil except + #wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation.
        String toBracketString(bool integer_mass, bool mass_delta, libcpp_vector[String] fixed_modifications) nogil except + # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation

        void setModification(Size index, const String& modification) nogil except + # wrap-doc:Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version

        void setModification(Size index, const ResidueModification& modification) nogil except + # wrap-doc:Sets the modification of AA at index by providing a ResidueModification object. Stricter than just looking for the name and adds the Modification to the DB if not present

        void setModificationByDiffMonoMass(Size index, double diffMonoMass) nogil except + # wrap-doc:Modifies the residue at index in the sequence and potentially in the ResidueDB

        void setNTerminalModification(String modification) nogil except + # wrap-doc:Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        void setNTerminalModification(const ResidueModification& mod) nogil except + # wrap-doc:Sets the N-terminal modification (copies and adds to database if not present)

        void setNTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term) nogil except + # wrap-doc:Sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)

        void setCTerminalModification(String modification) nogil except + # wrap-doc:Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        void setCTerminalModification(const ResidueModification& mod) nogil except + # wrap-doc:Sets the C-terminal modification (copies and adds to database if not present)

        void setCTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term) nogil except + # wrap-doc:Sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)

        String getNTerminalModificationName() nogil except + # wrap-doc:Returns the name (ID) of the N-terminal modification, or an empty string if none is set

        const ResidueModification * getNTerminalModification() nogil except + # wrap-doc:Returns a copy of the name N-terminal modification object, or None

        String getCTerminalModificationName() nogil except + # wrap-doc:Returns the name (ID) of the C-terminal modification, or an empty string if none is set

        const ResidueModification * getCTerminalModification() nogil except + # wrap-doc:Returns a copy of the name C-terminal modification object, or None

        # returns the residue at position index
        Residue getResidue(Size index) nogil except + # wrap-doc:Returns the residue at position index

        # returns the formula of the peptide
        EmpiricalFormula getFormula() nogil except + # wrap-doc:Convenience function with ResidueType=Full and charge = 0 by default
        EmpiricalFormula getFormula(ResidueType type_, Int charge) nogil except +

        # returns the average weight of the peptide
        double getAverageWeight() nogil except + # wrap-doc:Returns the average weight of the peptide
        double getAverageWeight(ResidueType type_, Int charge) nogil except +

        # returns the mono isotopic weight of the peptide
        double getMonoWeight() nogil except + # wrap-doc:Returns the mono isotopic weight of the peptide
        double getMonoWeight(ResidueType type_, Int charge) nogil except +

        # returns the mass-to-charge ratio of the peptide
        double getMZ(Int charge) nogil except + # wrap-doc:Returns the mass-to-charge ratio of the peptide
        double getMZ(Int charge, ResidueType type_) nogil except +

        # returns the number of residues
        Size size() nogil except + # wrap-doc:Returns the number of residues

        # returns a peptide sequence of the first index residues
        AASequence getPrefix(Size index) nogil except + # wrap-doc:Returns a peptide sequence of the first index residues

        # returns a peptide sequence of the last index residues
        AASequence getSuffix(Size index) nogil except + # wrap-doc:Returns a peptide sequence of the last index residues

        # returns a peptide sequence of number residues, beginning at position index
        AASequence getSubsequence(Size index, UInt number) nogil except + # wrap-doc:Returns a peptide sequence of number residues, beginning at position index

        # compute frequency table of amino acids
        void getAAFrequencies(Map[String, size_t]) nogil except + # wrap-ignore

        # returns true if the peptide contains the given residue
        bool has(Residue residue) nogil except + # wrap-doc:Returns true if the peptide contains the given residue

        # returns true if the peptide contains the given peptide
        # @note c-term and n-term mods are ignored
        bool hasSubsequence(AASequence peptide) nogil except + # wrap-doc:Returns true if the peptide contains the given peptide

        # returns true if the peptide has the given prefix
        # n-term mod is also checked (c-term as well, if prefix is of same length)
        bool hasPrefix(AASequence peptide) nogil except + # wrap-doc:Returns true if the peptide has the given prefix

        # returns true if the peptide has the given suffix
        # c-term mod is also checked (n-term as well, if suffix is of same length)
        bool hasSuffix(AASequence peptide) nogil except + # wrap-doc:Returns true if the peptide has the given suffix

        # predicate which is true if the peptide is N-term modified
        bool hasNTerminalModification() nogil except + # wrap-doc:Predicate which is true if the peptide is N-term modified

        # predicate which is true if the peptide is C-term modified
        bool hasCTerminalModification() nogil except + # wrap-doc:Predicate which is true if the peptide is C-term modified

        # returns true if any of the residues or termini are modified
        bool isModified() nogil except + # wrap-doc:Returns true if any of the residues or termini are modified

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS::AASequence":
        
        # static members
        AASequence fromString(String s, bool permissive) nogil except +  # wrap-attach:AASequence wrap-as:fromStringPermissive
        
        # static members
        AASequence fromString(String s) nogil except +  # wrap-attach:AASequence
