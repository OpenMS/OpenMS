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

#ifndef OPENMS_CHEMISTRY_RESIDUE_H
#define OPENMS_CHEMISTRY_RESIDUE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iosfwd>
#include <set>
#include <vector>

namespace OpenMS
{

  /**
      @ingroup Chemistry

      @brief Representation of a residue

      This class represents residues. Residues can have many different attributes, like
      the formula physico-chemical values of properties and so on.

      A very important property of residues are their modifications. By default no
      modification is present. Any modification which is present in the ModificationsDB can
      be applied, if appropriate.
  */
  class OPENMS_DLLAPI Residue
  {
public:

    /** @name Typedefs and Constants
    */
    //@{

    // Formulae that need to be added to the internal residues to get to fragment type
    // Formulae from http://www.matrixscience.com/help/fragmentation_help.html
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
      static const EmpiricalFormula to_full = getInternalToCTerm() + EmpiricalFormula("CO") - EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToYIon()
    {
      static const EmpiricalFormula to_full = getInternalToCTerm() + EmpiricalFormula("H");
      return to_full;
    }

    inline static const EmpiricalFormula& getInternalToZIon()
    {
      // Mind the "-"
      static const EmpiricalFormula to_full = getInternalToCTerm() - EmpiricalFormula("NH2");
      return to_full;
    }


    //@}

    /** @name Enums
    */
    //@{
    enum ResidueType
    {
      Full = 0,           // with N-terminus and C-terminus
      Internal,           // internal, without any termini
      NTerminal,           // only N-terminus
      CTerminal,           // only C-terminus
      AIon,           // N-terminus up to the C-alpha/carbonyl carbon bond
      BIon,           // N-terminus up to the peptide bond
      CIon,           // N-terminus up to the amide/C-alpha bond
      XIon,           // amide/C-alpha bond up to the C-terminus
      YIon,           // peptide bond up to the C-terminus
      ZIon,            // C-alpha/carbonyl carbon bond
      SizeOfResidueType
    };
    //@}

    /// returns the ion name given as a residue type
    static String getResidueTypeName(const ResidueType res_type);


    /** @name Constructors
    */
    //@{
    /// default constructor
    Residue();

    /// copy constructor
    Residue(const Residue & residue);

    /// detailed constructor
    Residue(const String & name,
            const String & three_letter_code,
            const String & one_letter_code,
            const EmpiricalFormula & formula);

    /// destructor
    virtual ~Residue();
    //@}

    /** @name Assignment
     */
    //@{
    /// assignment operator
    Residue & operator=(const Residue & residue);
    //@}

    /** Accessors
    */
    //@{
    /// sets the name of the residue
    void setName(const String & name);

    /// returns the name of the residue
    const String & getName() const;

    /// sets the short name of the residue, this name is used in the PeptideSequence for output
    void setShortName(const String & short_name);

    /// returns the short name of the residue
    const String & getShortName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<String> & synonyms);

    /// adds a synonym
    void addSynonym(const String & synonym);

    /// returns the synonyms
    const std::set<String> & getSynonyms() const;

    /// sets the name of the residue as three letter code
    void setThreeLetterCode(const String & three_letter_code);

    /// returns the name of the residue as three letter code
    const String & getThreeLetterCode() const;

    /// sets the name as one letter code
    void setOneLetterCode(const String & one_letter_code);

    /// returns the name as one letter code
    const String & getOneLetterCode() const;

    /// adds a neutral loss formula
    void addLossFormula(const EmpiricalFormula &);

    /// sets the neutral loss formulas
    void setLossFormulas(const std::vector<EmpiricalFormula> &);

    /// adds N-terminal losses
    void addNTermLossFormula(const EmpiricalFormula &);

    /// sets the N-terminal losses
    void setNTermLossFormulas(const std::vector<EmpiricalFormula> &);

    /// returns the neutral loss formulas
    const std::vector<EmpiricalFormula> & getLossFormulas() const;

    /// returns N-terminal loss formulas
    const std::vector<EmpiricalFormula> & getNTermLossFormulas() const;

    /// set the neutral loss molecule name
    void setLossNames(const std::vector<String> & name);

    /// sets the N-terminal loss names
    void setNTermLossNames(const std::vector<String> & name);

    /// add neutral loss molecule name
    void addLossName(const String & name);

    /// adds a N-terminal loss name
    void addNTermLossName(const String & name);

    /// gets neutral loss name (if there is one, else returns an empty string)
    const std::vector<String> & getLossNames() const;

    /// returns the N-terminal loss names
    const std::vector<String> & getNTermLossNames() const;

    /// set empirical formula of the residue (must be full, with N and C-terminus)
    void setFormula(const EmpiricalFormula & formula);

    /// returns the empirical formula of the residue
    EmpiricalFormula getFormula(ResidueType res_type = Full) const;

    /// sets average weight of the residue (must be full, with N and C-terminus)
    void setAverageWeight(double weight);

    /// returns average weight of the residue
    double getAverageWeight(ResidueType res_type = Full) const;

    /// sets mono weight of the residue (must be full, with N and C-terminus)
    void setMonoWeight(double weight);

    /// returns mono weight of the residue
    double getMonoWeight(ResidueType res_type = Full) const;

    /// sets by the name, this mod should be present in ModificationsDB
    void setModification(const String & name);

    /// returns the name of the modification to the modification
    const String & getModification() const;

    /// sets the low mass marker ions as a vector of formulas
    void setLowMassIons(const std::vector<EmpiricalFormula> & low_mass_ions);

    /// returns a vector of formulas with the low mass markers of the residue
    const std::vector<EmpiricalFormula> & getLowMassIons() const;

    /// sets the residue sets the amino acid is contained in
    void setResidueSets(const std::set<String> & residues_sets);

    /// adds a residue set to the residue sets
    void addResidueSet(const String & residue_sets);

    /// returns the residue sets this residue is contained in
    const std::set<String> & getResidueSets() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// true if the residue has neutral loss
    bool hasNeutralLoss() const;

    /// true if N-terminal neutral losses are set
    bool hasNTermNeutralLosses() const;

    /// equality operator
    bool operator==(const Residue & residue) const;

    /// inequality operator
    bool operator!=(const Residue & residue) const;

    /// equality operator for one letter code
    bool operator==(char one_letter_code) const;

    /// equality operator for one letter code
    bool operator!=(char one_letter_code) const;

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

    /// returns the backbone basicitiy if located in N-terminal direction
    double getBackboneBasicityLeft() const;

    /// sets the N-terminal direction backbone basicitiy
    void setBackboneBasicityLeft(double gb_bb_l);

    /// returns the C-terminal direction backbone basicitiy
    double getBackboneBasicityRight() const;

    /// sets the C-terminal direction backbone basicity
    void setBackboneBasicityRight(double gb_bb_r);

    /// true if the residue is a modified one
    bool isModified() const;

    /// true if the residue is contained in the set
    bool isInResidueSet(const String & residue_set);
    //@}

    /// ostream iterator to write the residue to a stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Residue & residue);

protected:

    // basic
    String name_;

    String short_name_;

    std::set<String> synonyms_;

    String three_letter_code_;

    String one_letter_code_;

    EmpiricalFormula formula_;

    EmpiricalFormula internal_formula_;

    double average_weight_;

    double mono_weight_;

    // modification
    bool is_modified_;

    String pre_mod_name_;

    String modification_;

    // loss
    std::vector<String> loss_names_;

    std::vector<EmpiricalFormula> loss_formulas_;

    std::vector<String> NTerm_loss_names_;

    std::vector<EmpiricalFormula> NTerm_loss_formulas_;

    double loss_average_weight_;

    double loss_mono_weight_;

    // low mass markers like immonium ions
    std::vector<EmpiricalFormula> low_mass_ions_;

    // pka values
    double pka_;

    // pkb values
    double pkb_;

    // pkc values
    double pkc_;

    double gb_sc_;

    double gb_bb_l_;

    double gb_bb_r_;

    // residue sets this amino acid is contained in
    std::set<String> residue_sets_;

  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Residue & residue);

}

#endif
