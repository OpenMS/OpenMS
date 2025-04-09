// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <set>

namespace OpenMS
{
  // forward declaration
  class Residue;

  /** @brief Representation of a modification on an amino acid residue

      This class represents a modification of a Residue. A residue modification
      has several attributes like the diff formula, a terminal specificity,
      a mass and maybe an origin which means a specific residue which it can
      be applied to. A residue modification can be represented by its Unimod name
      identifier (Id), e.g. "Oxidation (M)" or "Oxidation". This is a unique key which
      only occurs once in an OpenMS instance stored in the ModificationsDB.

      Example: methionine sulfoxide formation by oxidation of methionine

      @code
      Function              | Result
      ----------------------------------------------------
      getFullId()           | "Oxidation (M)"
      getId()               | "Oxidation"
      getFullName()         | "Oxidation or Hydroxylation"
      getUniModAccession()  | "UniMod:312"
      @endcode

      Note that some modifications are not explicitly defined from an input
      file but get added on the fly when reading amino acid sequences with
      bracket notation, e.g. "PEPTX[999]IDE". If there is no known modification
      corresponding to the indicated mass, then a new ResidueModification will
      be created which will return the initial string through "getFullId()" --
      which will either be "[999]" for internal modifications or ".[999]" for
      N/C-terminal modifications. Please use "isUserDefined" to check for
      user-defined modifications (those without 'Id' but with a 'FullId').

  */
  class OPENMS_DLLAPI ResidueModification
  {
public:

    /** Enums
    */
    //@{
    /** @brief Position where the modification is allowed to occur

            The allowed sites are
            Anywhere
            Any C-term
            Any N-term
            Protein C-term
            Protein N-term

            This does not describe which modifications are valid for a
            specific amino acid!
    */
    enum TermSpecificity
    {
      ANYWHERE = 0,        ///< The modification can go anywhere
      C_TERM = 1,          ///< The modification can only be C-terminal on the peptide
      N_TERM = 2,          ///< The modification can only be N-terminal on the peptide
      PROTEIN_C_TERM = 3,  ///< The modification can only be protein C-terminal (on a C-terminal peptide)
      PROTEIN_N_TERM = 4,  ///< The modification can only be protein N-terminal (on a N-terminal peptide)
      NUMBER_OF_TERM_SPECIFICITY
    };

    /** @brief Classification of the modification
    */
    enum SourceClassification
    {
      ARTIFACT = 0,           ///< A technical or chemical artefact
      HYPOTHETICAL,           ///< A hypothetical modification
      NATURAL,                ///< A natural modification
      POSTTRANSLATIONAL,      ///< A post-translational modification
      MULTIPLE,               ///< A combination of multiple modifications
      CHEMICAL_DERIVATIVE,    ///< A chemical derivative
      ISOTOPIC_LABEL,         ///< A chemical label (artificial)
      PRETRANSLATIONAL,       ///< A pre-translational modification
      OTHER_GLYCOSYLATION,    ///<  Other glycosylation (i.e. neither N nor O-linked)
      NLINKED_GLYCOSYLATION,  ///< N-linked glycosylation 
      AA_SUBSTITUTION,        ///< An amino acid substitution that presents like a modification
      OTHER,                  ///< Other modification
      NONSTANDARD_RESIDUE,    ///< A non-standard amino acid
      COTRANSLATIONAL,        ///< A co-translational modification
      OLINKED_GLYCOSYLATION,  ///< A O-linked glycosylation
      UNKNOWN,                ///< An unknown modification
      NUMBER_OF_SOURCE_CLASSIFICATIONS
    };
    //@}

    /** @name Constructors and Destructors
    */
    //@{

    /// Default constructor
    ResidueModification();

    /// Copy constructor
    ResidueModification(const ResidueModification&) = default;

    /// Move constructor
    ResidueModification(ResidueModification&&) = default;

    /// Destructor
    virtual ~ResidueModification();
    //@}

    /** @name Assignment operator
    */
    //@{

    /// Assignment operator
    ResidueModification& operator=(const ResidueModification&) = default;

    /// Move assignment operator
    ResidueModification& operator=(ResidueModification&&) & = default;
    //@}

    /** @name Accessors
    */
    //@{
    /// set the identifier of the modification
    void setId(const String& id);

    /// returns the identifier of the modification
    const String& getId() const;

    /**
       @brief Sets the full identifier (Unimod Accession + origin, if available)

       With empty argument, create a full ID based on (short) ID, terminal specificity and residue of origin.

       @throw Exception::MissingInformation if both argument @p full_id and member @p id_ are empty.
    */
    void setFullId(const String& full_id = "");

    /// returns the full identifier of the mod (Unimod accession + origin, if available)
    /**
       @brief Returns the full identifier (Unimod Accession + origin, if available)

       @note This field is used for user-defined modifications as well and will
       be set to a string such as "[999]" for internal modifications or
       ".[999]" for N/C-terminal modifications. Please use "isUserDefined" to
       check for user-defined modifications.

    */
    const String& getFullId() const;

    /// sets the unimod record id
    void setUniModRecordId(const Int& id);

    /// gets the unimod record id
    const Int& getUniModRecordId() const;

    /// returns the unimod accession if available
    const String getUniModAccession() const;

    /// set the MOD:XXXXX accession of PSI-MOD
    void setPSIMODAccession(const String& id);

    /// returns the PSI-MOD accession if available
    const String& getPSIMODAccession() const;

    /// sets the full name of the modification; must NOT contain the origin (or . for terminals!)
    void setFullName(const String& full_name);

    /// returns the full name of the modification
    const String& getFullName() const;

    /// sets the name of modification
    void setName(const String& name);

    /// returns the PSI-MS-label if available; e.g. Mascot uses this name
    const String& getName() const;

    /**
       @brief Sets the term specificity

       @throw Exception::InvalidValue if no valid specificity was given
    */
    void setTermSpecificity(TermSpecificity term_spec);

    /**
       @brief Sets the terminal specificity using a name

       Valid names: "C-term", "N-term", "none"

       @throw Exception::InvalidValue if no valid specificity was given
    */
    void setTermSpecificity(const String& name);

    /// returns terminal specificity
    TermSpecificity getTermSpecificity() const;

    /**
       @brief Returns the name of the terminal specificity

       By default, returns the name of the specificity set in member @p term_spec_.
       Alternatively, returns the name corresponding to argument @p term_spec.
    */
    String getTermSpecificityName(TermSpecificity term_spec = NUMBER_OF_TERM_SPECIFICITY) const;

    /**
       @brief Sets the origin (i.e. modified amino acid)

       @p origin must be a valid amino acid one-letter code or X, i.e. a letter from A to Y, excluding B and J.
       X represents any amino acid (for modifications with terminal specificity).

       @throw Exception::InvalidValue if @p origin is not in the valid range
    */
    void setOrigin(char origin);

    /// Returns the origin (i.e. modified amino acid)
    char getOrigin() const;

    /// classification as defined by the PSI-MOD
    void setSourceClassification(const String& classification);

    /// sets the source classification
    void setSourceClassification(SourceClassification classification);

    /// returns the source classification, if none was set, it is unspecific
    SourceClassification getSourceClassification() const;

    /// returns the classification
    String getSourceClassificationName(SourceClassification classification = NUMBER_OF_SOURCE_CLASSIFICATIONS) const;

    /// sets the average mass
    void setAverageMass(double mass);

    /// returns the average mass if set
    double getAverageMass() const;

    /// sets the monoisotopic mass (this must include the weight of the residue itself!)
    void setMonoMass(double mass);

    /// return the monoisotopic mass, or 0.0 if not set
    double getMonoMass() const;

    /// set the difference average mass
    void setDiffAverageMass(double mass);

    /// returns the difference average mass, or 0.0 if not set
    double getDiffAverageMass() const;

    /// sets the difference monoisotopic mass
    void setDiffMonoMass(double mass);

    /// returns the diff monoisotopic mass, or 0.0 if not set
    double getDiffMonoMass() const;

    /// set the formula (no masses will be changed)
    void setFormula(const String& composition);

    /// returns the chemical formula if set
    const String& getFormula() const;

    /// sets diff formula (no masses will be changed)
    void setDiffFormula(const EmpiricalFormula& diff_formula);

    /// returns the diff formula if one was set
    const EmpiricalFormula& getDiffFormula() const;

    /// sets the synonyms of that modification
    void setSynonyms(const std::set<String>& synonyms);

    /// adds a synonym to the unique list
    void addSynonym(const String& synonym);

    /// returns the set of synonyms
    const std::set<String>& getSynonyms() const;

    /// sets the neutral loss formula
    void setNeutralLossDiffFormulas(const std::vector<EmpiricalFormula>& diff_formulas);

    /// returns the neutral loss diff formula (if available)
    const std::vector<EmpiricalFormula>& getNeutralLossDiffFormulas() const;

    /// set the neutral loss mono weight
    void setNeutralLossMonoMasses(std::vector<double> mono_masses);

    /// returns the neutral loss mono weight
    std::vector<double> getNeutralLossMonoMasses() const;

    /// set the neutral loss average weight
    void setNeutralLossAverageMasses(std::vector<double> average_masses);

    /// returns the neutral loss average weight
    std::vector<double> getNeutralLossAverageMasses() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if a neutral loss formula is set
    bool hasNeutralLoss() const;

    /// returns true if it is a user-defined modification (empty id)
    bool isUserDefined() const;

    /// equality operator
    bool operator==(const ResidueModification& modification) const;

    /// inequality operator
    bool operator!=(const ResidueModification& modification) const;

    /// less operator
    bool operator<(const ResidueModification& modification) const;

    //@}

    /// Creates a new modification from a mass and adds it to ModificationsDB.
    /// If not terminal, needs a Residue to be put on.

    /// @param mod The mass to put between the brackets (might contain +/- at the front)
    /// @param mass Basically, the same as mod, just as double (since usually both representations are present when calling this function and to avoid overhead??)
    /// @param delta_mass Is the given mass a delta mass (i.e. does @p mod contain a = or -)?
    /// @param specificity To which site can this mod be applied?
    /// @param residue [only required for ANYWHERE term spec] Residue with further information (e.g. residue weights) for the new mod
    /// @return a new or existing mod; registered to ModDB in both cases, so the pointer is non-owning (FullId is e.g. M[+1234.1] and FullName [+1234.1]. Id and Name are empty as defined for a "user-defined" mod.
    static const ResidueModification* createUnknownFromMassString(const String& mod,
                                                                  const double mass,
                                                                  const bool delta_mass,
                                                                  const TermSpecificity specificity,
                                                                  const Residue* residue = nullptr);

    /** @brief Merge a set of mods to a given modification (usually the one which is already present, but can be null)

    If only one mod is combined in total, it is not changed to an unknown mod but remains a 'known' mod.
    If base is already contained in @p addons, it is not added again.
    
    All mods given here must have the same term specificity and origin (which might be 'X', i.e. no restriction), otherwise a Precondition exception is thrown.

    If base and addons is empty, a null_ptr is returned.

    @param base An already present mod, can be a nullptr
    @param addons A set of mods to add on top of the mod. 
    @param allow_unknown_masses If any input (incl. base) is already an unknown mass, nothing is done
    @param residue [only required for ANYWHERE term spec] Residue with further information (e.g. residue weights) for the new mod
    @return A (new custom) mod, which is registered in ModificationsDB if needed.
    @throws Exception::Precondition if term spec or origins to not match between all given mods
    **/
    static const ResidueModification* combineMods(const ResidueModification* base,
                                                  const std::set<const ResidueModification*>& addons,
                                                  bool allow_unknown_masses = false,
                                                  const Residue* residue = nullptr);

    /// Convert to string (incl. origin/terminal), in order of preference:
    ///  + using the ID, as X(ID), e.g. 'M(Oxidation)' or '.(Acetyl)'
    ///  + using the FullName
    ///  + using the delta_mono_mass, e.g. 'M[+15.65]'
    ///  + using the mono_mass, e.g. 'M[56.23]'
    ///
    /// The mono_mass must not be negative (undistinguishable to delta_mono_mass when parsing)
    String toString() const;

    /// converts the mass to a string with preceding '+' or '-' sign
    /// e.g. '-19.34' or '+1.003'
    static String getDiffMonoMassString(const double diff_mono_mass);

    /// return a string of the form '[+&gt;mass&lt;] (the '+' might be a '-', if mass is negative).
    static String getDiffMonoMassWithBracket(const double diff_mono_mass);

    /// return a string of the form '[&gt;mass&lt;]
    static String getMonoMassWithBracket(const double mono_mass);

protected:
    String id_;

    String full_id_;

    String psi_mod_accession_;

    // The UniMod record id (an integer)
    Int unimod_record_id_;

    String full_name_;

    String name_;

    TermSpecificity term_spec_;

    char origin_;

    SourceClassification classification_;

    double average_mass_;

    double mono_mass_;

    double diff_average_mass_;

    double diff_mono_mass_;

    String formula_;

    EmpiricalFormula diff_formula_;

    std::set<String> synonyms_;

    std::vector<EmpiricalFormula> neutral_loss_diff_formulas_;

    std::vector<double> neutral_loss_mono_masses_;

    std::vector<double> neutral_loss_average_masses_;
  };
}
