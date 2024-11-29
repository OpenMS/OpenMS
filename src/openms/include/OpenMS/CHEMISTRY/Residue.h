// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Jang Jang Jin$
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>



#include <array>
#include <iosfwd>
#include <set>
#include <vector>

namespace OpenMS
{
  class ResidueModification;

  /**
      @ingroup Chemistry

      @brief Representation of an amino acid residue.

      This class represents an amino acid. These can have many different attributes, like
      the formula physico-chemical values of properties and so on.

      A very important property of amino acid residues are their modifications. By default no
      modification is present. Any modification which is present in the ModificationsDB can
      be applied, if appropriate.
  */
  class OPENMS_DLLAPI Residue
  {
    friend class ResidueDB;

public:

    /** @name Formula conversion
     *
     * @brief Computes empirical formula required to add to the desired type
     *
     * Computes the empirical formula required to be added to convert an
     * internal residue (inside an AA sequence) to that of a residue of the
     * desired type. For example, to obtain the conversion formula for an
     * internal ion to a "y ion", use getInternalToYTerm().
     *
     * Formulae that need to be added to the internal residues to get to
     * fragment type from http://www.matrixscience.com/help/fragmentation_help.html
     * 
     * Which means that we follow Biemann nomenclature.
     * For a small description of the differernt ion types, see the enum @ref ResidueType.
     */
    //@{

    inline static const EmpiricalFormula& getInternalToFull()
    {
      static const EmpiricalFormula to_full = EmpiricalFormula("H2O");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToNTerm()
    {
      static const EmpiricalFormula to_full = EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToCTerm()
    {
      static const EmpiricalFormula to_full = EmpiricalFormula("OH");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToAIon()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToNTerm() - EmpiricalFormula("CHO");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToBIon()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToNTerm() - EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToCIon()
    {
      static const EmpiricalFormula to_full =
        getInternalToNTerm() + EmpiricalFormula("NH2");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToXIon()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToCTerm() + EmpiricalFormula("CO") - EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToYIon()
    {
      static const EmpiricalFormula to_full =
        getInternalToCTerm() + EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToZIon()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToCTerm() - EmpiricalFormula("NH2");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToZp1Ion()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToCTerm() - EmpiricalFormula("NH");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToZp2Ion()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full =
        getInternalToCTerm() - EmpiricalFormula("N");
      return to_full;
    }

    //@}

    /** @name Enums
    */
    //@{

    /// Residue types: Note that all weights and elemental compositions of fragment "ions" are given for their neutral forms.
    /// Furthermore, all fragment ion types are based on the Biemann nomenclature (http://www.matrixscience.com/help/fragmentation_help.html)
    /// See https://github.com/OpenMS/OpenMS/issues/7219 for a discussion with more details and links.
    enum ResidueType
    {
      Full = 0,       ///< with N-terminus and C-terminus
      Internal,       ///< internal residue, without any termini
      NTerminal,      ///< only N-terminus
      CTerminal,      ///< only C-terminus
      AIon,           ///< MS:1001229 N-terminus up to the C-alpha/carbonyl carbon bond
      BIon,           ///< MS:1001224 N-terminus up to the peptide bond
      CIon,           ///< MS:1001231 N-terminus up to the amide/C-alpha bond
      XIon,           ///< MS:1001228 amide/C-alpha bond up to the C-terminus
      YIon,           ///< MS:1001220 peptide bond up to the C-terminus
      ZIon,           ///< MS:1001230 C-alpha/carbonyl carbon bond [CID fragment]
      Zp1Ion,         ///< MS:1001230 C-alpha/carbonyl carbon bond (free radical, z+1 "ion") [main EAD fragment]
      Zp2Ion,         ///< MS:1001230 C-alpha/carbonyl carbon bond (free radical, z+2 "ion" with additional abstracted hydrogen) [EAD fragment at higher precursor charges]
      Precursor,      ///< MS:1001523 Precursor ion
      BIonMinusH20,   ///< MS:1001222 b ion without water
      YIonMinusH20,   ///< MS:1001223 y ion without water
      BIonMinusNH3,   ///< MS:1001232 b ion without ammonia
      YIonMinusNH3,   ///< MS:1001233 y ion without ammonia
      NonIdentified,  ///< MS:1001240 Non-identified ion
      Unannotated,    ///< no stored annotation
      SizeOfResidueType
    };
    //@}
    
    /// Names corresponding to the ResidueType enum
    static inline std::array<std::string_view, Residue::ResidueType::SizeOfResidueType> names_of_residuetype {
      "full",
      "internal",
      "N-terminal",
      "C-terminal",
      "a-ion",
      "b-ion",
      "c-ion",
      "x-ion",
      "y-ion",
      "z-ion",
      "z+1-ion",
      "z+2-ion",
      "precursor-ion",
      "b-H2O-ion",
      "y-H2O-ion",
      "b-NH3-ion",
      "y-NH3-ion",
      "Non-identified ion",
      "unannotated"
    };

    /// returns the ion name given as a residue type
    static String getResidueTypeName(const ResidueType res_type);

    /** @name Constructors
    */
    //@{

    /// Default constructor (needed by pyOpenMS)
    Residue();

    /// Copy constructor
    Residue(const Residue&) = default;

    /// Move constructor
    Residue(Residue&&) = default;
           
    // Detailed constructor 
    Residue(const String& name,
            const String& three_letter_code,
            const String& one_letter_code,
            const EmpiricalFormula& formula,
            double pka = 0,
            double pkb = 0,
            double pkc = -1,
            double gb_sc = 0,
            double gb_bb_l = 0,
            double gb_bb_r = 0,
            const std::set<String>& synonyms = std::set<String>());

    /// Destructor
    virtual ~Residue();
    //@}

    /** @name Assignment
     */
    //@{

    /// Assignment operator
    Residue& operator=(const Residue&) = default;

    /// Move assignment operator
    Residue& operator=(Residue&&) & = default;
    //@}

    /** @name Accessors
    */
    //@{
    /// sets the name of the residue
    void setName(const String& name);

    /// returns the name of the residue
    const String& getName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<String>& synonyms);

    /// adds a synonym
    void addSynonym(const String& synonym);

    /// returns the synonyms
    const std::set<String>& getSynonyms() const;

    /// sets the name of the residue as three letter code (String of size 3)
    void setThreeLetterCode(const String& three_letter_code);

    /// returns the name of the residue as three letter code (String of size 3)
    const String& getThreeLetterCode() const;

    /// sets the name as one letter code (String of size 1)
    void setOneLetterCode(const String& one_letter_code);

    /// returns the name as one letter code (String of size 1)
    const String& getOneLetterCode() const;

    /// adds a neutral loss formula
    void addLossFormula(const EmpiricalFormula&);

    /// sets the neutral loss formulas
    void setLossFormulas(const std::vector<EmpiricalFormula>&);

    /// adds N-terminal losses
    void addNTermLossFormula(const EmpiricalFormula&);

    /// sets the N-terminal losses
    void setNTermLossFormulas(const std::vector<EmpiricalFormula>&);

    /// returns the neutral loss formulas
    const std::vector<EmpiricalFormula>& getLossFormulas() const;

    /// returns N-terminal loss formulas
    const std::vector<EmpiricalFormula>& getNTermLossFormulas() const;

    /// set the neutral loss molecule name
    void setLossNames(const std::vector<String>& name);

    /// sets the N-terminal loss names
    void setNTermLossNames(const std::vector<String>& name);

    /// add neutral loss molecule name
    void addLossName(const String& name);

    /// adds a N-terminal loss name
    void addNTermLossName(const String& name);

    /// gets neutral loss name (if there is one, else returns an empty string)
    const std::vector<String>& getLossNames() const;

    /// returns the N-terminal loss names
    const std::vector<String>& getNTermLossNames() const;

    /// set empirical formula of the residue (must be full, with N and C-terminus)
    void setFormula(const EmpiricalFormula& formula);

    /// returns the empirical formula of the residue
    EmpiricalFormula getFormula(ResidueType res_type = Full) const;

    /// sets average weight of the residue (must be full, with N and C-terminus)
    void setAverageWeight(double weight);

    /// returns average weight of the residue
    double getAverageWeight(ResidueType res_type = Full) const;

    /// sets monoisotopic weight of the residue (must be full, with N and C-terminus)
    void setMonoWeight(double weight);

    /// returns monoisotopic weight of the residue
    double getMonoWeight(ResidueType res_type = Full) const;

    /// returns a pointer to the modification, or a null pointer if none is set
    const ResidueModification* getModification() const;

    /// sets the modification by name; the mod should be present in ModificationsDB
    void setModification(const String& name);

    /// sets the modification by existing ResMod (make sure it exists in ModificationsDB)
    void setModification(const ResidueModification* mod);

    /// sets the modification by looking for an exact match in the DB first, otherwise creating a
    /// new entry
    void setModification(const ResidueModification& mod);

    /// sets a modification by monoisotopic mass difference. Searches in DBs first with a tolerance.
    /// If not found, creates a new entry with name = OneLetterResidueCode[+/-diffMonoMass] and adds this as user-defined mod.
    void setModificationByDiffMonoMass(double diffMonoMass);
    
    /// returns the name (ID) of the modification, or an empty string if none is set
    const String& getModificationName() const;

    /// sets the low mass marker ions as a vector of formulas
    void setLowMassIons(const std::vector<EmpiricalFormula>& low_mass_ions);

    /// returns a vector of formulas with the low mass markers of the residue
    const std::vector<EmpiricalFormula>& getLowMassIons() const;

    /// sets the residue sets the amino acid is contained in (e.g. Natural20)
    void setResidueSets(const std::set<String>& residues_sets);

    /// adds a residue set to the residue sets (e.g. Natural20)
    void addResidueSet(const String& residue_sets);

    /// returns the residue sets this residue is contained in (e.g. Natural20)
    const std::set<String>& getResidueSets() const;

    /// returns the pka of the residue
    double getPka() const;

    /// returns the pkb of the residue
    double getPkb() const;

    /// returns the pkc of the residue if it exists otherwise -1
    double getPkc() const;

    /// calculates the isoelectric point using the pk* values
    double getPiValue() const;

    /// sets the pka of the residue
    void setPka(double value);

    /// sets the pkb of the residue
    void setPkb(double value);

    /// sets the pkc of the residue
    void setPkc(double value);

    /// returns the side chain basicity
    double getSideChainBasicity() const;

    /// sets the side chain basicity
    void setSideChainBasicity(double gb_sc);

    /// returns the backbone basicity if located in N-terminal direction
    double getBackboneBasicityLeft() const;

    /// sets the N-terminal direction backbone basicity
    void setBackboneBasicityLeft(double gb_bb_l);

    /// returns the C-terminal direction backbone basicity
    double getBackboneBasicityRight() const;

    /// sets the C-terminal direction backbone basicity
    void setBackboneBasicityRight(double gb_bb_r);
    //@}

    /** @name Predicates
    */
    //@{
    /// true if the residue has neutral loss
    bool hasNeutralLoss() const;

    /// true if N-terminal neutral losses are set
    bool hasNTermNeutralLosses() const;

    /// equality operator
    bool operator==(const Residue& residue) const;

    /// inequality operator
    bool operator!=(const Residue& residue) const;

    /// equality operator for one letter code
    bool operator==(char one_letter_code) const;

    /// equality operator for one letter code
    bool operator!=(char one_letter_code) const;

    /// true if the residue is a modified one
    bool isModified() const;

    /// true if the residue is contained in the set
    bool isInResidueSet(const String& residue_set);
    //@}

    /// helper for mapping residue types to letters for Text annotations and labels
    static std::string residueTypeToIonLetter(const ResidueType& res_type);

    /// Write as Origin+Modification, e.g. M(Oxidation), or X[945.34] or N[+14.54] for user-defined mods.
    /// This requires the Residue to have a valid OneLetterCode and an optional (but valid) ResidueModification (see ResidueModification::toString())
    String toString() const;

    /// ostream iterator to write the residue to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Residue& residue);

protected:

    /// the name of the residue
    String name_ = "unknown";

    std::set<String> synonyms_;

    String three_letter_code_;

    String one_letter_code_;

    EmpiricalFormula formula_;

    EmpiricalFormula internal_formula_;

    double average_weight_ = 0;

    double mono_weight_ = 0;

    /// pointer to the modification 
    const ResidueModification* modification_ = nullptr;

    // loss
    std::vector<String> loss_names_;

    std::vector<EmpiricalFormula> loss_formulas_;

    std::vector<String> NTerm_loss_names_;

    std::vector<EmpiricalFormula> NTerm_loss_formulas_;

    double loss_average_weight_;

    double loss_mono_weight_;

    /// low mass markers like immonium ions
    std::vector<EmpiricalFormula> low_mass_ions_;

    // pka values
    double pka_ = 0;

    // pkb values
    double pkb_ = 0;

    // pkc values
    double pkc_ = -1.0;

    /// SideChainBasicity
    double gb_sc_ = 0;

    /// BackboneBasicityLeft
    double gb_bb_l_ = 0;

    /// BackboneBasicityRight
    double gb_bb_r_ = 0;

    /// residue sets this amino acid is contained in
    std::set<String> residue_sets_;

    // pre-calculated residue type delta weights for more efficient weight calculation
    static const double internal_to_full_monoweight_;
    static const double internal_to_nterm_monoweight_;
    static const double internal_to_cterm_monoweight_;
    static const double internal_to_a_monoweight_;
    static const double internal_to_b_monoweight_;
    static const double internal_to_c_monoweight_;
    static const double internal_to_x_monoweight_;
    static const double internal_to_y_monoweight_;
    static const double internal_to_z_monoweight_;
    static const double internal_to_zp1_monoweight_;
    static const double internal_to_zp2_monoweight_;
  };

  // write 'name threelettercode onelettercode formula'
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Residue& residue);

}
