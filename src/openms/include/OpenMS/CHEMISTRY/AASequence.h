// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <vector>
#include <iosfwd>
#include <map>

namespace OpenMS
{

  /**
      @brief Representation of a peptide/protein sequence

      This class represents amino acid sequences in %OpenMS. An AASequence
      instance primarily contains a sequence of residues. The sequence is
      represented as a vector of pointers to instances of Residue. Each amino
      acid has only one instance, which is accessible using the ResidueDB
      instance (singleton).

      To create an AASequence instance for a specific amino acid sequence, use
      the AASequence::fromString function. For example,
      <tt>AASequence::fromString(".DFPIANGER.")</tt> produces an instance of
      AASequence for the peptide "DFPIANGER". Please note that both the N- and
      the C-terminal are explicitly represented by dots.

      A critical property of amino acid sequences is that they can be modified.
      Which means that one or more amino acids are chemically modified, e.g.
      oxidized. This is represented via Residue instances which carry a
      ResidueModification object. This is also handled in the ResidueDB.

      Modifications are specified using a unique string identifier present in
      the ModificationsDB in round brackets after the modified amino acid or by
      providing the mass of the residue in square brackets. For example
      <tt>AASequence::fromString(".DFPIAM(Oxidation)GER.")</tt> creates an
      instance of the peptide "DFPIAMGER" with an oxidized methionine
      (<tt>AASequence::fromString(".DFPIAM(UniMod:35)GER.")</tt>,
      <tt>AASequence::fromString(".DFPIAM[+16]GER.")</tt> and
      <tt>AASequence::fromString(".DFPIAM[147]GER.")</tt> are all equivalent).
      N- and C-terminal modifications are represented by brackets to the right
      of the dots terminating the sequence. For example,
      <tt>".(Dimethyl)DFPIAMGER."</tt> and <tt>".DFPIAMGER.(Label:18O(2))"</tt>
      represent the labelling of the N- and C-terminus respectively, but
      <tt>".DFPIAMGER(Phospho)."</tt> will be interpreted as a phosphorylation
      of the last arginine at its side chain.

      Note there is a subtle difference between
      <tt>AASequence::fromString(".DFPIAM[+16]GER.")</tt> and
      <tt>AASequence::fromString(".DFPIAM[+15.9949]GER.")</tt> -- while the
      former will try to find the @e first modification matching to a mass
      difference of 16 +/- 0.5, the latter will try to find the @e closest
      matching modification to the exact mass. This usually gives the intended
      results while the first approach may not.

      Arbitrary/unknown amino acids (usually due to an unknown modification)
      can be specified using tags preceded by X: "X[weight]". This indicates a
      new amino acid ("X") with the specified weight, e.g. "RX[148.5]T"". Note
      that this tag does not alter the amino acids to the left (R) or right
      (T).  Rather, X represents an amino acid on its own. Be careful when
      converting such AASequence objects to an EmpiricalFormula using
      getFormula(), as tags will not be considered in this case (there exists
      no formula for them).  However, they have an influence on getMonoWeight()
      and getAverageWeight()!

      @note For C/N terminal modifications, the absolute mass is assumed to be
      1 (H) for the N-terminus and 17 (OH) for the C-terminus, therefore a
      modification specified as absolute n[43]PEPTIDE would translate to
      n[+42]PEPTIDE using relative masses. Note that there can be ambiguity in
      cases where the loss includes the terminal amino acids.

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI AASequence final
  {
public:

    class Iterator;

    /** @brief ConstIterator for AASequence

        AASequence constant iterator
    */
    class OPENMS_DLLAPI ConstIterator final
    {
    public:
      // TODO Iterator constructor for ConstIterator

      typedef const Residue& const_reference;
      typedef Residue& reference;
      typedef const Residue* const_pointer;
      typedef std::vector<const Residue*>::difference_type difference_type;
      typedef Residue value_type;
      typedef const Residue* pointer;
      typedef std::random_access_iterator_tag iterator_category;

      /** @name Constructors and destructors
       */
      //@{
      /// default constructor
      ConstIterator() = default;

      /// detailed constructor with pointer to the vector and offset position
      ConstIterator(const std::vector<const Residue*>* vec_ptr, difference_type position)
        : vector_{vec_ptr},
          position_{position}
      {
      }

      /// copy constructor
      ConstIterator(const ConstIterator& rhs) = default;

      /// copy constructor from Iterator
      ConstIterator(const AASequence::Iterator& rhs) :
        vector_(rhs.vector_),
        position_(rhs.position_)
      {
      }

      /// destructor
      ~ConstIterator() = default;

      //@}

      /// assignment operator
      ConstIterator& operator=(const ConstIterator& rhs) = default;

      /** @name Operators
      */
      //@{
      /// dereference operator
      const_reference operator*() const
      {
        return *(*vector_)[position_];
      }

      /// dereference operator
      const_pointer operator->() const
      {
        return (*vector_)[position_];
      }

      /// forward jump operator
      const ConstIterator operator+(difference_type diff) const
      {
        return ConstIterator(vector_, position_ + diff);
      }

      difference_type operator-(ConstIterator rhs) const
      {
        return position_ - rhs.position_;
      }

      /// backward jump operator
      const ConstIterator operator-(difference_type diff) const
      {
        return ConstIterator(vector_, position_ - diff);
      }

      /// equality comparator
      bool operator==(const ConstIterator& rhs) const
      {
        return vector_ == rhs.vector_ && position_ == rhs.position_;
      }

      /// inequality operator
      bool operator!=(const ConstIterator& rhs) const
      {
        return vector_ != rhs.vector_ || position_ != rhs.position_;
      }

      /// increment operator
      ConstIterator& operator++()
      {
        ++position_;
        return *this;
      }

      /// decrement operator
      ConstIterator& operator--()
      {
        --position_;
        return *this;
      }

      //@}

protected:

      // pointer to the AASequence vector
      const std::vector<const Residue*>* vector_ {};

      // position in the AASequence vector
      difference_type position_ {};
    };


    /** @brief Iterator class for AASequence

        Mutable iterator for AASequence
    */
    class OPENMS_DLLAPI Iterator final
    {
public:

      friend class AASequence::ConstIterator;

      typedef const Residue& const_reference;
      typedef Residue& reference;
      typedef const Residue* const_pointer;
      typedef const Residue* pointer;
      typedef std::vector<const Residue*>::difference_type difference_type;

      /** @name Constructors and destructors
      */
      //@{
      /// default constructor
      Iterator() = default;

      /// detailed constructor with pointer to the vector and offset position
      Iterator(std::vector<const Residue*>* vec_ptr, difference_type position)
        : vector_ {vec_ptr},
          position_{position}
      {
      }

      /// copy constructor
      Iterator(const Iterator& rhs) = default;

      /// destructor
      ~Iterator() = default;

      //@}

      /// assignment operator
      Iterator& operator=(const Iterator& rhs)
      {
        if (this != &rhs)
        {
          position_ = rhs.position_;
          vector_ = rhs.vector_;
        }
        return *this;
      }

      /** @name Operators
      */
      //@{
      /// dereference operator
      const_reference operator*() const
      {
        return *(*vector_)[position_];
      }

      /// dereference operator
      const_pointer operator->() const
      {
        return (*vector_)[position_];
      }

      /// mutable dereference operator
      pointer operator->()
      {
        return (*vector_)[position_];
      }

      /// forward jump operator
      const Iterator operator+(difference_type diff) const
      {
        return Iterator(vector_, position_ + diff);
      }

      difference_type operator-(Iterator rhs) const
      {
        return position_ - rhs.position_;
      }

      /// backward jump operator
      const Iterator operator-(difference_type diff) const
      {
        return Iterator(vector_, position_ - diff);
      }

      /// equality comparator
      bool operator==(const Iterator& rhs) const
      {
        return vector_ == rhs.vector_ && position_ == rhs.position_;
      }

      /// inequality operator
      bool operator!=(const Iterator& rhs) const
      {
        return vector_ != rhs.vector_ || position_ != rhs.position_;
      }

      /// increment operator
      Iterator& operator++()
      {
        ++position_;
        return *this;
      }

      /// decrement operator
      Iterator& operator--()
      {
        --position_;
        return *this;
      }

      //@}

protected:

      // pointer to the AASequence vector
      std::vector<const Residue*>* vector_ {};

      // position in the AASequence vector
      difference_type position_ {};
    };

    /** @name Constructors and Destructors
    */
    //@{

    /// Default constructor
    AASequence() = default;

    /// Copy constructor
    AASequence(const AASequence&) = default;

    /// Move constructor
    AASequence(AASequence&&) = default;

    /// Destructor
    ~AASequence() = default;
    //@}

    /// Assignment operator
    AASequence& operator=(const AASequence&) = default;

    /// Move assignment operator
    AASequence& operator=(AASequence&&) = default;

    /// check if sequence is empty
    bool empty() const;

    /** @name Accessors
    */
    //@{

    /**
        @brief returns the peptide as string with modifications embedded in brackets

        Uses round brackets when possible (id is known) or square brackets for
        unknown modifications where only the mass is known.

        i.e.: .[43]PEPC(Carbamidomethyl)PEPM[147]PEPR.[-1]

        @note For unknown modifications, the function will attempt to use the
        exact same format used in the input
    */
    String toString() const;

    /// returns the peptide as string without any modifications or (e.g., "PEPTIDER")
    String toUnmodifiedString() const;

    /**
        @brief returns the peptide as string with UniMod-style modifications embedded in brackets

        Annotates modification with UniMod identifier (when identifier is
        known) and uses square brackets for unknown modifications (only mass is known).

        i.e.: .[43]PEPC(UniMod:4)PEPM[147]PEPR.[16]
    */
    String toUniModString() const;

    /**
        @brief create a TPP compatible string of the modified sequence using bracket notation.

        Instead of using the modification names, it writes the modification masses in brackets

        i.e.:

         - n[43]PEPC[160]PEPM[147]PEPRc[16]
         - n[+42]PEPC[+57]PEPM[+16]PEPRc[-1]

        will be produced, depending on whether relative or absolute masses are used.

        @param integer_mass Whether to use integer masses in brackets (default is true, if false, accurate masses will be written)
        @param mass_delta Whether to write absolute masses M[147] or relative mass deltas M[+16] (default is false)
        @param fixed_modifications Optional list of fixed modifications that should not be added to the output (they are considered to be present in all cases)

        @note Using integer masses may mean that there could be multiple modifications mapping to the same mass
    */
    String toBracketString(bool integer_mass = true,
                           bool mass_delta = false,
                           const std::vector<String> & fixed_modifications = std::vector<String>()) const;

    /// set the modification of the residue at position index.
    /// if an empty string is passed replaces the residue with its unmodified version
    void setModification(Size index, const String& modification);

    /// sets the modification of AA at @p index by providing an already, potentially modified residue
    void setModification(Size index, const Residue* modification);

    /// sets the modification of AA at @p index by providing a pointer to a ResidueModification object found in the ModificationsDB
    void setModification(Size index, const ResidueModification* modification);

    /// sets the modification of AA at @p index by providing a ResidueModification object
    /// stricter than just looking for the name and adds the Modification to the DB if not present
    void setModification(Size index, const ResidueModification& modification);

    /// modifies the residue at @p index in the sequence and potentially in the ResidueDB
    void setModificationByDiffMonoMass(Size index, double diffMonoMass);

    /// sets the N-terminal modification (by lookup in the mod names of the ModificationsDB)
    /// throws if nothing is found (since the name is not enough information to create a new mod)
    void setNTerminalModification(const String& modification);

    /// sets the N-terminal modification
    void setNTerminalModification(const ResidueModification* modification);

    /// sets the N-terminal modification (copies and adds to database if not present)
    void setNTerminalModification(const ResidueModification& mod);

    /// sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
    void setNTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term);

    /// returns the name (ID) of the N-terminal modification, or an empty string if none is set
    const String& getNTerminalModificationName() const;

    /// returns a pointer to the N-terminal modification, or zero if none is set
    const ResidueModification* getNTerminalModification() const;

    /// sets the C-terminal modification (by lookup in the mod names of the ModificationsDB)
    /// throws if nothing is found (since the name is not enough information to create a new mod)
    void setCTerminalModification(const String& modification);

    /// sets the C-terminal modification (must be present in the database)
    void setCTerminalModification(const ResidueModification* modification);

    /// sets the C-terminal modification (copies and adds to database if not present)
    void setCTerminalModification(const ResidueModification& mod);

    /// sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
    void setCTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term);

    /// returns the name (ID) of the C-terminal modification, or an empty string if none is set
    const String& getCTerminalModificationName() const;

    /// returns a pointer to the C-terminal modification, or zero if none is set
    const ResidueModification* getCTerminalModification() const;

    /// returns a pointer to the residue at position @p index
    const Residue& getResidue(Size index) const;

    /// returns the formula of the peptide
    EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns the average weight of the peptide
    double getAverageWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns the mono isotopic weight of the peptide in the given ionic form
    /// @note will not (and cannot) control whether the required ion can exist
    /// (e.g. x/c ions for monomers) as it does not do fragmentation but rather
    /// supplementing/deduction of the sequence to its ionic form.
    double getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

    /// returns mass-to-charge ratio of the peptide in the given ionic form
    /// @note will not (and cannot) control whether the required ion can exist
    /// (e.g. x/c ions for monomers) as it does not do fragmentation but rather
    /// supplementing/deduction of the sequence to its ionic form.
    /// @throws Exception::InvalidValue if @p charge==0
    double getMZ(Int charge, Residue::ResidueType type = Residue::Full) const;

    /// returns a pointer to the residue at given position
    const Residue& operator[](Size index) const;

    /// adds the residues of the peptide
    AASequence operator+(const AASequence& peptide) const;

    /// adds the residues of a peptide
    AASequence& operator+=(const AASequence&);

    /// adds the residues of the peptide
    AASequence operator+(const Residue* residue) const;

    /// adds the residues of a peptide
    AASequence& operator+=(const Residue*);

    /// returns the number of residues
    Size size() const;

    /// returns a peptide sequence of the first index residues
    AASequence getPrefix(Size index) const;

    /// returns a peptide sequence of the last index residues
    AASequence getSuffix(Size index) const;

    /// returns a peptide sequence of number residues, beginning at position index
    AASequence getSubsequence(Size index, UInt number) const;

    /// compute frequency table of amino acids
    void getAAFrequencies(std::map<String, Size>& frequency_table) const;

    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the peptide contains the given residue
    bool has(const Residue& residue) const;

    /// returns true if the peptide contains the given peptide
    /// @note c-term and n-term mods are ignored
    bool hasSubsequence(const AASequence& peptide) const;

    /// returns true if the peptide has the given prefix
    /// n-term mod is also checked (c-term as well, if prefix is of same length)
    bool hasPrefix(const AASequence& peptide) const;

    /// returns true if the peptide has the given suffix
    /// c-term mod is also checked (n-term as well, if suffix is of same length)
    bool hasSuffix(const AASequence& peptide) const;

    /// predicate which is true if the peptide is N-term modified
    bool hasNTerminalModification() const;

    /// predicate which is true if the peptide is C-term modified
    bool hasCTerminalModification() const;

    /// returns true if any of the residues or termini are modified
    bool isModified() const;

    /// equality operator. Two sequences are equal iff all amino acids including PTMs are equal
    bool operator==(const AASequence& rhs) const;

    /// lesser than operator which compares the C-term mods, sequence including PTMS and N-term mods; can be used for maps
    bool operator<(const AASequence& rhs) const;

    /// inequality operator. Complement of equality operator.
    bool operator!=(const AASequence& rhs) const;
    //@}

    /** @name Iterators
    */
    //@{
    inline Iterator begin() { return Iterator(&peptide_, 0); }

    inline ConstIterator begin() const { return ConstIterator(&peptide_, 0); }

    inline Iterator end() { return Iterator(&peptide_, (Int) peptide_.size()); }

    inline ConstIterator end() const { return ConstIterator(&peptide_, (Int) peptide_.size()); }
    //@}

    /** @name Stream operators
    */
    //@{
    /// writes a peptide to an output stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AASequence& peptide);

    /// reads a peptide from an input stream
    friend OPENMS_DLLAPI std::istream& operator>>(std::istream& is, const AASequence& peptide);
    //@}

    /**
      @brief create AASequence object by parsing an OpenMS string

      @param s Input string
      @param permissive If set, skip spaces and replace stop codon symbols ("*", "#", "+") by "X" (unknown amino acid) during parsing

      @throws Exception::ParseError if an invalid string representation of an AA sequence is passed
    */
    static AASequence fromString(const String& s,
                                 bool permissive = true);

    /**
      @brief create AASequence object by parsing a C string (character array)

      @param s Input string
      @param permissive If set, skip spaces and replace stop codon symbols ("*", "#", "+") by "X" (unknown amino acid) during parsing

      @throws Exception::ParseError if an invalid string representation of an AA sequence is passed
    */
    static AASequence fromString(const char* s,
                                 bool permissive = true);

    /// @brief constructor from String
    /// @param s A String representing the amino acid sequence
    explicit AASequence(const String& s);

    /// @brief constructor from C string
    /// @param s A C-style string representing the amino acid sequence
    explicit AASequence(const char* s);

    /// @brief constructor from String
    /// @param s A String representing the amino acid sequence
    /// @param permissive If set, skip spaces and replace stop codon symbols ("*", "#", "+") by "X" (unknown amino acid) during parsing
    explicit AASequence(const String& s, bool permissive);

    /// @brief constructor from C string
    /// @param s A C-style string representing the amino acid sequence
    /// @param permissive If set, skip spaces and replace stop codon symbols ("*", "#", "+") by "X" (unknown amino acid) during parsing
    explicit AASequence(const char* s, bool permissive);

  protected:

    std::vector<const Residue*> peptide_;

    const ResidueModification* n_term_mod_ = nullptr;

    const ResidueModification* c_term_mod_ = nullptr;

    /**
      @brief Parses modifications in round brackets (an identifier)

      If dot notation is used it resolves cterm ambiguity based on the presence
      of the dot.

      @param str_it Current position in the string to be parsed
      @param str Full input string
      @param aas Current AASequence object (will be modified with the correct residue added)
      @param specificity Whether the current modification should be interpreted as N- or C-terminal

      @return Position at which to continue parsing
    */
    static String::ConstIterator parseModRoundBrackets_(const String::ConstIterator str_it,
                                                        const String& str,
                                                        AASequence& aas,
                                                        const ResidueModification::TermSpecificity& specificity);

    /**
      @brief Parses modifications in square brackets (a mass)

      If dot notation is used it resolves cterm ambiguity based on the presence
      of the dot.

      @param str_it Current position in the string to be parsed
      @param str Full input string
      @param aas Current AASequence object (will be modified with the correct residue added)
      @param specificity Whether the current modification should be interpreted as N- or C-terminal

      @return Position at which to continue parsing
    */
    static String::ConstIterator parseModSquareBrackets_(const String::ConstIterator str_it,
                                                         const String& str,
                                                         AASequence& aas,
                                                         const ResidueModification::TermSpecificity& specificity);

    static void parseString_(const String& peptide, AASequence& aas,
                             bool permissive = true);
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AASequence& peptide);

  OPENMS_DLLAPI std::istream& operator>>(std::istream& os, const AASequence& peptide);

} // namespace OpenMS


