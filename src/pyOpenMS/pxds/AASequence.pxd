from libcpp cimport bool
from libcpp.map cimport map as libcpp_map
from Types cimport *
from String cimport *
from Residue cimport *
from EmpiricalFormula cimport *
from ResidueModification cimport *

cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS":

    cdef cppclass AASequence:
        # wrap-hash:
        #  toString().c_str()
        #
        # wrap-doc:
        #  Representation of a peptide/protein sequence
        #  This class represents amino acid sequences in OpenMS. An AASequence
        #  instance primarily contains a sequence of residues. 

        AASequence() except + nogil 
        AASequence(AASequence &) except + nogil
    
        AASequence(const String&) except + nogil # wrap-doc:Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE")
        AASequence(const String&, bool permissive) except + nogil # wrap-doc:Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE"), permissive allows for '+', '*', and '#' in the sequence

        AASequence operator+(AASequence) except + nogil 
        AASequence iadd(AASequence) except + nogil  # wrap-as:operator+=

        bool operator==(AASequence &rhs) except + nogil
        bool operator!=(AASequence &rhs) except + nogil

        # Note that this is a const-ref, so we cannot easily set residues
        Residue operator[](size_t) except + nogil  # wrap-upper-limit:size()

        # check if sequence is empty
        bool empty() except + nogil  # wrap-doc:Check if sequence is empty

        # returns the peptide as string with modifications embedded in brackets
        String toString() except + nogil  # wrap-doc:Returns the peptide as string with modifications embedded in brackets

        # returns the peptide as string without any modifications
        String toUnmodifiedString() except + nogil  # wrap-doc:Returns the peptide as string without any modifications

        # returns the peptide as string without any modifications
        String toUniModString() except + nogil  # wrap-doc:Returns the peptide as string with UniMod-style modifications embedded in brackets

        String toBracketString() except + nogil  # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation. Uses integer mass by default
        String toBracketString(bool integer_mass) except + nogil  # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation
        String toBracketString(bool integer_mass, bool mass_delta) except + nogil  #wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation.
        String toBracketString(bool integer_mass, bool mass_delta, libcpp_vector[String] fixed_modifications) except + nogil  # wrap-doc:Create a TPP compatible string of the modified sequence using bracket notation

        void setModification(Size index, const String& modification) except + nogil  # wrap-doc:Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version

        void setModification(Size index, const ResidueModification& modification) except + nogil  # wrap-doc:Sets the modification of AA at index by providing a ResidueModification object. Stricter than just looking for the name and adds the Modification to the DB if not present

        void setModificationByDiffMonoMass(Size index, double diffMonoMass) except + nogil  # wrap-doc:Modifies the residue at index in the sequence and potentially in the ResidueDB

        void setNTerminalModification(String modification) except + nogil  # wrap-doc:Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        void setNTerminalModification(const ResidueModification& mod) except + nogil  # wrap-doc:Sets the N-terminal modification (copies and adds to database if not present)

        void setNTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term) except + nogil  # wrap-doc:Sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)

        void setCTerminalModification(String modification) except + nogil  # wrap-doc:Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        void setCTerminalModification(const ResidueModification& mod) except + nogil  # wrap-doc:Sets the C-terminal modification (copies and adds to database if not present)

        void setCTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term) except + nogil  # wrap-doc:Sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)

        String getNTerminalModificationName() except + nogil  # wrap-doc:Returns the name (ID) of the N-terminal modification, or an empty string if none is set

        const ResidueModification * getNTerminalModification() except + nogil  # wrap-doc:Returns a copy of the name N-terminal modification object, or None

        String getCTerminalModificationName() except + nogil  # wrap-doc:Returns the name (ID) of the C-terminal modification, or an empty string if none is set

        const ResidueModification * getCTerminalModification() except + nogil  # wrap-doc:Returns a copy of the name C-terminal modification object, or None

        # returns the residue at position index
        Residue getResidue(Size index) except + nogil  # wrap-doc:Returns the residue at position index

        # returns the formula of the peptide
        EmpiricalFormula getFormula() except + nogil  # wrap-doc:Convenience function with ResidueType=Full and charge = 0 by default
        EmpiricalFormula getFormula(ResidueType type_, Int charge) except + nogil 

        # returns the average weight of the peptide
        double getAverageWeight() except + nogil  # wrap-doc:Returns the average weight of the peptide
        double getAverageWeight(ResidueType type_, Int charge) except + nogil 

        # returns the mono isotopic weight of the peptide
        double getMonoWeight() except + nogil  # wrap-doc:Returns the mono isotopic weight of the peptide
        double getMonoWeight(ResidueType type_, Int charge) except + nogil 

        # returns the mass-to-charge ratio of the peptide
        double getMZ(Int charge) except + nogil  # wrap-doc:Returns the mass-to-charge ratio of the peptide
        double getMZ(Int charge, ResidueType type_) except + nogil 

        # returns the number of residues
        Size size() except + nogil  # wrap-doc:Returns the number of residues

        # returns a peptide sequence of the first index residues
        AASequence getPrefix(Size index) except + nogil  # wrap-doc:Returns a peptide sequence of the first index residues

        # returns a peptide sequence of the last index residues
        AASequence getSuffix(Size index) except + nogil  # wrap-doc:Returns a peptide sequence of the last index residues

        # returns a peptide sequence of number residues, beginning at position index
        AASequence getSubsequence(Size index, UInt number) except + nogil  # wrap-doc:Returns a peptide sequence of number residues, beginning at position index

        # compute frequency table of amino acids
        void getAAFrequencies(libcpp_map[String, size_t]) except + nogil  # wrap-ignore

        # returns true if the peptide contains the given residue
        bool has(Residue residue) except + nogil  # wrap-doc:Returns true if the peptide contains the given residue

        # returns true if the peptide contains the given peptide
        # @note c-term and n-term mods are ignored
        bool hasSubsequence(AASequence peptide) except + nogil  # wrap-doc:Returns true if the peptide contains the given peptide

        # returns true if the peptide has the given prefix
        # n-term mod is also checked (c-term as well, if prefix is of same length)
        bool hasPrefix(AASequence peptide) except + nogil  # wrap-doc:Returns true if the peptide has the given prefix

        # returns true if the peptide has the given suffix
        # c-term mod is also checked (n-term as well, if suffix is of same length)
        bool hasSuffix(AASequence peptide) except + nogil  # wrap-doc:Returns true if the peptide has the given suffix

        # predicate which is true if the peptide is N-term modified
        bool hasNTerminalModification() except + nogil  # wrap-doc:Predicate which is true if the peptide is N-term modified

        # predicate which is true if the peptide is C-term modified
        bool hasCTerminalModification() except + nogil  # wrap-doc:Predicate which is true if the peptide is C-term modified

        # returns true if any of the residues or termini are modified
        bool isModified() except + nogil  # wrap-doc:Returns true if any of the residues or termini are modified

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/AASequence.h>" namespace "OpenMS::AASequence":
        
        # static members
        AASequence fromString(String s, bool permissive) except + nogil   # wrap-attach:AASequence wrap-as:fromStringPermissive wrap-doc:deprecated. Use AASequence(String, bool) instead.
        
        # static members
        AASequence fromString(String s) except + nogil   # wrap-attach:AASequence wrap-doc:deprecated. Use AASequence(String) instead.
