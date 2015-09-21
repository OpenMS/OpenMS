// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H
#define OPENMS_CHEMISTRY_RESIDUEMODIFICATION_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <set>

namespace OpenMS
{
  // forward declaration
  class Residue;

  /** @brief Representation of a modification

          This class represents a modification of a residue. A residue modification
          has several attributes like the diff formula, a terminal specificity
          a mass and maybe an origin which means a specific residue which it can
          be applied to. A residue modification can be represented by its UniMod name
          identifier, e.g. "Oxidation (M)" or "Oxidation". This is a unique key which
          only occurs once in an OpenMS instance stored in the ModificationsDB.

          Example: methionine sulfoxide formation by oxidation of methionine

          getFullId()           = "Oxidation (M)"
          getId()               = "Oxidation"
          getFullName()         = "Oxidation or Hydroxylation"
          getUniModAccession()  = "UniMod:312"
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

            This does not describe the amino acids which are valid for a
            specific amino acid!

    */
    enum Term_Specificity
    {
      ANYWHERE = 0,
      C_TERM = 1,
      N_TERM = 2,
      PROTEIN_C_TERM = 3,
      PROTEIN_N_TERM = 4,
      NUMBER_OF_TERM_SPECIFICITY
    };

    /** @brief Classification of the modification

    */
    enum Source_Classification
    {
      ARTIFACT = 0,
      HYPOTHETICAL,
      NATURAL,
      POSTTRANSLATIONAL,
      MULTIPLE,
      CHEMICAL_DERIVATIVE,
      ISOTOPIC_LABEL,
      PRETRANSLATIONAL,
      OTHER_GLYCOSYLATION,
      NLINKED_GLYCOSYLATION,
      AA_SUBSTITUTION,
      OTHER,
      NONSTANDARD_RESIDUE,
      COTRANSLATIONAL,
      OLINKED_GLYCOSYLATION,
      UNKNOWN,
      NUMBER_OF_SOURCE_CLASSIFICATIONS
    };
    //@}



    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    ResidueModification();

    /// copy constructor
    ResidueModification(const ResidueModification & modification);

    /// destructor
    virtual ~ResidueModification();
    //@}

    /** @name Assignment operator
    */
    //@{
    /// assignment operator
    ResidueModification & operator=(const ResidueModification & modification);
    //@}

    /** @name Accessors
    */
    //@{
    /// set the identifier of the modification
    void setId(const String & id);

    /// returns the identifier of the modification
    const String & getId() const;

    /// set full identifier (UniMod Accession + origin, if available)
    void setFullId(const String & full_id);

    /// returns the full id of the mod (UniMod accession + origin, if available)
    const String & getFullId() const;

    /// sets the unimod accession
    void setUniModAccession(const String & id);

    /// returns the unimod accession if available
    const String & getUniModAccession() const;

    /// set the MOD:XXXXX accession of PSI-MOD
    void setPSIMODAccession(const String & id);

    /// returns the PSI-MOD accession if available
    const String & getPSIMODAccession() const;

    /// sets the full name of the modification
    void setFullName(const String & full_name);

    /// returns the full name of the modification
    const String & getFullName() const;

    /// sets the name of modification
    void setName(const String & name);

    /// returns the PSI-MS-label if available; e.g. Mascot uses this name
    const String & getName() const;

    /// sets the term specificity
    void setTermSpecificity(Term_Specificity term_spec);

    /// sets the terminal specificity using a name (valid: "C-term","N-term","none")
    void setTermSpecificity(const String & name);

    /// returns terminal specificity
    Term_Specificity getTermSpecificity() const;

    /// returns the terminal specificity name which is set or given as parameter
    String getTermSpecificityName(Term_Specificity = NUMBER_OF_TERM_SPECIFICITY) const;

    ///sets the origin (i.e. amino acid)
    void setOrigin(const String & origin);

    /// returns the origin (i.e. amino acid) if set
    const String & getOrigin() const;

    /// classification as defined by the PSI-MOD
    void setSourceClassification(const String & classification);

    /// sets the source classification
    void setSourceClassification(Source_Classification classification);

    /// returns the source classification, if none was set, it is unspecific
    Source_Classification getSourceClassification() const;

    /// returns the classification
    String getSourceClassificationName(Source_Classification classification = NUMBER_OF_SOURCE_CLASSIFICATIONS) const;

    /// sets the average mass
    void setAverageMass(double mass);

    /// returns the average mass if set
    double getAverageMass() const;

    /// sets the monoisotopic mass
    void setMonoMass(double mass);

    /// return the monoisotopic mass, if set
    double getMonoMass() const;

    /// set the difference average mass
    void setDiffAverageMass(double mass);

    /// returns the difference average mass if set
    double getDiffAverageMass() const;

    /// sets the difference monoisotopic mass
    void setDiffMonoMass(double mass);

    /// returns the diff monoisotopic mass if set
    double getDiffMonoMass() const;

    /// set the formula
    void setFormula(const String & composition);

    /// returns the chemical formula if set
    const String & getFormula() const;

    /// sets diff formula
    void setDiffFormula(const EmpiricalFormula & diff_formula);

    /// returns the diff formula if one was set
    const EmpiricalFormula & getDiffFormula() const;

    /// sets the synonyms of that modification
    void setSynonyms(const std::set<String> & synonyms);

    /// adds a synonym to the unique list
    void addSynonym(const String & synonym);

    /// returns the set of synonyms
    const std::set<String> & getSynonyms() const;

    /// sets the neutral loss formula
    void setNeutralLossDiffFormula(const EmpiricalFormula & loss);

    /// returns the neutral loss diff formula (if available)
    const EmpiricalFormula & getNeutralLossDiffFormula() const;

    /// set the neutral loss mono weight
    void setNeutralLossMonoMass(double mono_mass);

    /// returns the neutral loss mono weight
    double getNeutralLossMonoMass() const;

    /// set the neutral loss average weight
    void setNeutralLossAverageMass(double average_mass);

    /// returns the neutral loss average weight
    double getNeutralLossAverageMass() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if a neutral loss formula is set
    bool hasNeutralLoss() const;

    /// equality operator
    bool operator==(const ResidueModification & modification) const;

    /// inequality operator
    bool operator!=(const ResidueModification & modification) const;
    //@}

protected:

    String id_;

    String full_id_;

    String psi_mod_accession_;

    String unimod_accession_;

    String full_name_;

    String name_;

    Term_Specificity term_spec_;

    String origin_;

    Source_Classification classification_;

    double average_mass_;

    double mono_mass_;

    double diff_average_mass_;

    double diff_mono_mass_;

    String formula_;

    EmpiricalFormula diff_formula_;

    std::set<String> synonyms_;

    EmpiricalFormula neutral_loss_diff_formula_;

    double neutral_loss_mono_mass_;

    double neutral_loss_average_mass_;
  };
}

#endif
