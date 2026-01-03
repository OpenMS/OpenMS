// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <functional>
#include <iosfwd>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Representation of a ribonucleotide (modified or unmodified)

      The available information is based on the Modomics database (http://modomics.genesilico.pl/modifications/).

      @see RibonucleotideDB
  */
  class OPENMS_DLLAPI Ribonucleotide
  {
    friend class RibonucleotideDB;
  public:

    enum TermSpecificityNuc
    {
      ANYWHERE,
      FIVE_PRIME,
      THREE_PRIME,
      NUMBER_OF_TERM_SPECIFICITY
    };

    /** @name Constructors
    */
    //@{
    /// Constructor
    Ribonucleotide(const String& name = "unknown ribonucleotide",
                   const String& code = ".",
                   const String& new_code = "",
                   const String& html_code = ".",
                   const EmpiricalFormula& formula = EmpiricalFormula(),
                   char origin = '.',
                   double mono_mass = 0.0,
                   double avg_mass = 0.0,
                   enum TermSpecificityNuc term_spec = ANYWHERE,
                   const EmpiricalFormula& baseloss_formula =
                   default_baseloss_);

    /// Copy constructor
    Ribonucleotide(const Ribonucleotide& ribo) = default;

    /// Destructor
    virtual ~Ribonucleotide();
    //@}

    /** @name Assignment
     */
    //@{
    /// assignment operator
    Ribonucleotide& operator=(const Ribonucleotide& ribo) = default;
    //@}


    /** @name Equality
     */
    //@{
    /// Equality operator
    bool operator==(const Ribonucleotide& ribonucleotide) const;
    //@}

    /** Accessors
     */
    //@{

    /// Return the short name
    const String getCode() const;

    /// Set the short name
    void setCode(const String& code);

    /// Get the name of the ribonucleotide
    const String getName() const;

    /// Set the name of the ribonucleotide
    void setName(const String& name);

    /// Get formula for the ribonucleotide
    const EmpiricalFormula getFormula() const;

    /// Set the empirical formula for the ribonucleotide
    void setFormula(const EmpiricalFormula& formula);

    /// Get the monoisotopic mass of the ribonucleotide
    double getMonoMass() const;

    /// Set the monoisotopic mass of the ribonucleotide
    void setMonoMass(double mono_mass);

    /// Set the average mass of the ribonucleotide
    double getAvgMass() const;

    /// Get the average mass of the ribonucleotide
    void setAvgMass(double avg_mass);

    /// Get the "new" (Modomics) code
    const String getNewCode() const;

    /// Set the "new" (Modomics) code
    void setNewCode(const String &new_code);

    /// ostream iterator to write the residue to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Ribonucleotide& ribo);

    /// Get the code of the unmodified base (e.g., "A", "C", ...)
    char getOrigin() const;

    /// Set the code of the unmodified base (e.g., "A", "C", ...)
    void setOrigin(char origin);

    /// Set the HTML (RNAMods) code
    String getHTMLCode() const;

    /// Get the HTML (RNAMods) code
    void setHTMLCode(const String& html_code);

    /// Get the terminal specificity
    enum TermSpecificityNuc getTermSpecificity() const;

    /// Set the terminal specificity
    void setTermSpecificity(enum TermSpecificityNuc term_spec);

    /// Get sum formula after loss of the nucleobase
    const EmpiricalFormula getBaselossFormula() const;

    /// Set the sum formula after loss of the nucleobase
    void setBaselossFormula(const EmpiricalFormula& formula);

    //@}

    /// Return true if this is a modified ribonucleotide and false otherwise
    bool isModified() const;

    /// Return whether this is an "ambiguous" modification (representing isobaric modifications on the base/ribose)
    bool isAmbiguous() const;

  protected:
    /// Default value for sum formula after nucleobase loss
    static const EmpiricalFormula default_baseloss_;

    String name_; ///< full name
    String code_; ///< short name
    String new_code_; ///< Modomics code
    String html_code_; ///< RNAMods code
    EmpiricalFormula formula_; ///< sum formula
    char origin_; ///< character of unmodified version of ribonucleotide
    double mono_mass_; ///< monoisotopic mass
    double avg_mass_; ///< average mass
    enum TermSpecificityNuc term_spec_; ///< terminal specificity
    EmpiricalFormula baseloss_formula_; ///< sum formula after loss of the nucleobase
  };

  /// Stream output operator
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Ribonucleotide& ribo);

  /// Dummy nucleotide used to represent 5' and 3' chain ends. Usually, just the phosphates.
  using RibonucleotideChainEnd = Ribonucleotide;

} // namespace OpenMS

// Hash function specialization for Ribonucleotide
// Placed in std namespace to allow use with std::unordered_map/set
namespace std
{
  /**
   * @brief Hash function for OpenMS::Ribonucleotide.
   *
   * Computes a hash based on all fields used in operator==:
   * name_, code_, new_code_, html_code_, formula_, origin_,
   * mono_mass_, avg_mass_, term_spec_, and baseloss_formula_.
   *
   * @note Hash is consistent with operator==: equal objects produce equal hashes.
   */
  template<>
  struct hash<OpenMS::Ribonucleotide>
  {
    std::size_t operator()(const OpenMS::Ribonucleotide& ribo) const noexcept
    {
      std::size_t seed = 0;

      // Hash all fields used in operator==
      // String fields: name_, code_, new_code_, html_code_
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ribo.getName()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ribo.getCode()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ribo.getNewCode()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(ribo.getHTMLCode()));

      // EmpiricalFormula fields: formula_, baseloss_formula_
      OpenMS::hash_combine(seed, std::hash<OpenMS::EmpiricalFormula>{}(ribo.getFormula()));
      OpenMS::hash_combine(seed, std::hash<OpenMS::EmpiricalFormula>{}(ribo.getBaselossFormula()));

      // Primitive fields: origin_ (char), mono_mass_, avg_mass_ (double), term_spec_ (enum)
      OpenMS::hash_combine(seed, OpenMS::hash_char(ribo.getOrigin()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(ribo.getMonoMass()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(ribo.getAvgMass()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(ribo.getTermSpecificity())));

      return seed;
    }
  };
} // namespace std
