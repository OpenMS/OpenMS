// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/CHEMISTRY/ResidueModification.h>

using namespace std;

namespace OpenMS
{

  ResidueModification::ResidueModification() :
    term_spec_(ResidueModification::ANYWHERE),
    classification_(ResidueModification::ARTIFACT),
    average_mass_(0.0),
    mono_mass_(0.0),
    diff_average_mass_(0.0),
    diff_mono_mass_(0.0),
    neutral_loss_mono_mass_(0.0),
    neutral_loss_average_mass_(0.0)
  {
  }

  ResidueModification::ResidueModification(const ResidueModification & rhs) :
    id_(rhs.id_),
    full_id_(rhs.full_id_),
    psi_mod_accession_(rhs.psi_mod_accession_),
    unimod_accession_(rhs.unimod_accession_),
    full_name_(rhs.full_name_),
    name_(rhs.name_),
    term_spec_(rhs.term_spec_),
    origin_(rhs.origin_),
    classification_(rhs.classification_),
    average_mass_(rhs.average_mass_),
    mono_mass_(rhs.mono_mass_),
    diff_average_mass_(rhs.diff_average_mass_),
    diff_mono_mass_(rhs.diff_mono_mass_),
    formula_(rhs.formula_),
    diff_formula_(rhs.diff_formula_),
    synonyms_(rhs.synonyms_),
    neutral_loss_diff_formula_(rhs.neutral_loss_diff_formula_),
    neutral_loss_mono_mass_(rhs.neutral_loss_mono_mass_),
    neutral_loss_average_mass_(rhs.neutral_loss_average_mass_)
  {
  }

  ResidueModification & ResidueModification::operator=(const ResidueModification & rhs)
  {
    if (this != &rhs)
    {
      id_ = rhs.id_;
      full_id_ = rhs.full_id_;
      psi_mod_accession_ = rhs.psi_mod_accession_;
      unimod_accession_ = rhs.unimod_accession_;
      full_name_ = rhs.full_name_;
      name_ = rhs.name_;
      term_spec_ = rhs.term_spec_;
      origin_ = rhs.origin_;
      classification_ = rhs.classification_;
      average_mass_ = rhs.average_mass_;
      mono_mass_ = rhs.mono_mass_;
      diff_average_mass_ = rhs.diff_average_mass_;
      diff_mono_mass_ = rhs.diff_mono_mass_;
      formula_ = rhs.formula_;
      diff_formula_ = rhs.diff_formula_;
      synonyms_ = rhs.synonyms_;
      neutral_loss_diff_formula_ = rhs.neutral_loss_diff_formula_;
      neutral_loss_mono_mass_ = rhs.neutral_loss_mono_mass_;
      neutral_loss_average_mass_ = rhs.neutral_loss_average_mass_;
    }

    return *this;
  }

  bool ResidueModification::operator==(const ResidueModification & rhs) const
  {
    return id_ == rhs.id_ &&
           full_id_ == rhs.full_id_ &&
           psi_mod_accession_ == rhs.psi_mod_accession_ &&
           unimod_accession_ == rhs.unimod_accession_ &&
           full_name_ == rhs.full_name_ &&
           name_ == rhs.name_ &&
           term_spec_ == rhs.term_spec_ &&
           origin_ == rhs.origin_ &&
           classification_ == rhs.classification_ &&
           average_mass_ == rhs.average_mass_ &&
           mono_mass_ == rhs.mono_mass_ &&
           diff_average_mass_ == rhs.diff_average_mass_ &&
           diff_mono_mass_ == rhs.diff_mono_mass_ &&
           formula_ == rhs.formula_ &&
           diff_formula_ == rhs.diff_formula_ &&
           synonyms_ == rhs.synonyms_ &&
           neutral_loss_diff_formula_ == rhs.neutral_loss_diff_formula_ &&
           neutral_loss_mono_mass_ == rhs.neutral_loss_mono_mass_ &&
           neutral_loss_average_mass_ == rhs.neutral_loss_average_mass_;
  }

  bool ResidueModification::operator!=(const ResidueModification & rhs) const
  {
    return !(*this == rhs);
  }

  ResidueModification::~ResidueModification()
  {

  }

  void ResidueModification::setId(const String & id)
  {
    id_ = id;
  }

  const String & ResidueModification::getId() const
  {
    return id_;
  }

  void ResidueModification::setFullId(const String & full_id)
  {
    full_id_ = full_id;
  }

  const String & ResidueModification::getFullId() const
  {
    return full_id_;
  }

  void ResidueModification::setPSIMODAccession(const String & id)
  {
    psi_mod_accession_ = id;
  }

  const String & ResidueModification::getPSIMODAccession() const
  {
    return psi_mod_accession_;
  }

  void ResidueModification::setUniModAccession(const String & id)
  {
    unimod_accession_ = id;
  }

  const String & ResidueModification::getUniModAccession() const
  {
    return unimod_accession_;
  }

  void ResidueModification::setFullName(const String & full_name)
  {
    full_name_ = full_name;
  }

  const String & ResidueModification::getFullName() const
  {
    return full_name_;
  }

  void ResidueModification::setName(const String & name)
  {
    name_ = name;
  }

  const String & ResidueModification::getName() const
  {
    return name_;
  }

  void ResidueModification::setTermSpecificity(Term_Specificity term_spec)
  {
    term_spec_ = term_spec;
  }

  void ResidueModification::setTermSpecificity(const String & term_spec)
  {
    if (term_spec == "C-term")
    {
      term_spec_ = C_TERM;
      return;
    }
    if (term_spec == "N-term")
    {
      term_spec_ = N_TERM;
      return;
    }
    if (term_spec == "none")
    {
      term_spec_ = ANYWHERE;
      return;
    }
    cerr << "ResidueModification: cannot convert '" << term_spec << "' into term specificity!" << endl;
    return;
  }

  ResidueModification::Term_Specificity ResidueModification::getTermSpecificity() const
  {
    return term_spec_;
  }

  String ResidueModification::getTermSpecificityName(Term_Specificity term_spec) const
  {
    if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
    {
      term_spec = term_spec_;
    }
    switch (term_spec)
    {
    case C_TERM: return "C-term";

    case N_TERM: return "N-term";

    default:         // ANYWHERE
      if (term_spec != ANYWHERE)
      {
        cerr << "ResidueModification: cannot convert '" << term_spec << "' into term specificity name!" << endl;
      }
      return "none";
    }
  }

  void ResidueModification::setOrigin(const String & origin)
  {
    origin_ = origin;
  }

  const String & ResidueModification::getOrigin() const
  {
    return origin_;
  }

  void ResidueModification::setSourceClassification(Source_Classification classification)
  {
    classification_ = classification;
  }

  void ResidueModification::setSourceClassification(const String & classification)
  {
    String c = classification;
    c.toLower();
    if (c == "artifact" || c == "artefact") // unimod uses Artefact (BE) not Artifcat (AE)
    {
      classification_ = ARTIFACT;
      return;
    }
    if (c == "natural")
    {
      classification_ = NATURAL;
      return;
    }
    if (c == "hypothetical")
    {
      classification_ = HYPOTHETICAL;
      return;
    }
    if (c == "post-translational")
    {
      classification_ = POSTTRANSLATIONAL;
      return;
    }
    if (c == "multiple")
    {
      classification_ = MULTIPLE;
      return;
    }
    if (c == "chemical derivative")
    {
      classification_ = CHEMICAL_DERIVATIVE;
      return;
    }
    if (c == "isotopic label")
    {
      classification_ = ISOTOPIC_LABEL;
      return;
    }
    if (c == "pre-translational")
    {
      classification_ = PRETRANSLATIONAL;
      return;
    }
    if (c == "other glycosylation")
    {
      classification_ = OTHER_GLYCOSYLATION;
      return;
    }
    if (c == "n-linked glycosylation")
    {
      classification_ = NLINKED_GLYCOSYLATION;
      return;
    }
    if (c == "aa substitution")
    {
      classification_ = AA_SUBSTITUTION;
      return;
    }
    if (c == "other")
    {
      classification_ = OTHER;
      return;
    }
    if (c == "non-standard residue")
    {
      classification_ = NONSTANDARD_RESIDUE;
      return;
    }
    if (c == "co-translational")
    {
      classification_ = COTRANSLATIONAL;
      return;
    }
    if (c == "o-linked glycosylation")
    {
      classification_ = OLINKED_GLYCOSYLATION;
      return;
    }

    classification_ = UNKNOWN;

    //cerr << "ResidueModification: Unknown source classification '" << classification << "'" << endl;
    return;
  }

  ResidueModification::Source_Classification ResidueModification::getSourceClassification() const
  {
    return classification_;
  }

  String ResidueModification::getSourceClassificationName(Source_Classification classification) const
  {
    if (classification == NUMBER_OF_SOURCE_CLASSIFICATIONS)
    {
      classification = classification_;
    }
    switch (classification)
    {
    case ARTIFACT: return "Artefact"; // return Artefact (BE) not Artifcat (AE)

    case NATURAL:  return "Natural";

    case HYPOTHETICAL: return "Hypothetical";

    case POSTTRANSLATIONAL: return "Post-translational";

    case MULTIPLE:  return "Multiple";

    case CHEMICAL_DERIVATIVE: return "Chemical derivative";

    case ISOTOPIC_LABEL: return "Isotopic label";

    case PRETRANSLATIONAL:  return "Pre-translational";

    case OTHER_GLYCOSYLATION: return "Other glycosylation";

    case NLINKED_GLYCOSYLATION: return "N-linked glycosylation";

    case AA_SUBSTITUTION:  return "AA substitution";

    case OTHER: return "Other";

    case NONSTANDARD_RESIDUE: return "Non-standard residue";

    case COTRANSLATIONAL: return "Co-translational";

    case OLINKED_GLYCOSYLATION:  return "O-linked glycosylation";

    case UNKNOWN:  return "";

    default: return "Unknown";
    }
    return "Unknown";
  }

  void ResidueModification::setAverageMass(DoubleReal mass)
  {
    average_mass_ = mass;
  }

  DoubleReal ResidueModification::getAverageMass() const
  {
    return average_mass_;
  }

  void ResidueModification::setMonoMass(DoubleReal mass)
  {
    mono_mass_ = mass;
  }

  DoubleReal ResidueModification::getMonoMass() const
  {
    return mono_mass_;
  }

  void ResidueModification::setDiffAverageMass(DoubleReal mass)
  {
    diff_average_mass_ = mass;
  }

  DoubleReal ResidueModification::getDiffAverageMass() const
  {
    return diff_average_mass_;
  }

  void ResidueModification::setDiffMonoMass(DoubleReal mass)
  {
    diff_mono_mass_ = mass;
  }

  DoubleReal ResidueModification::getDiffMonoMass() const
  {
    return diff_mono_mass_;
  }

  void ResidueModification::setFormula(const String & formula)
  {
    formula_ = formula;
  }

  const String & ResidueModification::getFormula() const
  {
    return formula_;
  }

  void ResidueModification::setDiffFormula(const EmpiricalFormula & diff_formula)
  {
    diff_formula_ = diff_formula;
  }

  const EmpiricalFormula & ResidueModification::getDiffFormula() const
  {
    return diff_formula_;
  }

  void ResidueModification::addSynonym(const String & synonym)
  {
    synonyms_.insert(synonym);
  }

  void ResidueModification::setSynonyms(const set<String> & synonyms)
  {
    synonyms_ = synonyms;
  }

  const set<String> & ResidueModification::getSynonyms() const
  {
    return synonyms_;
  }

  void ResidueModification::setNeutralLossDiffFormula(const EmpiricalFormula & diff_formula)
  {
    neutral_loss_diff_formula_ = diff_formula;
  }

  const EmpiricalFormula & ResidueModification::getNeutralLossDiffFormula() const
  {
    return neutral_loss_diff_formula_;
  }

  void ResidueModification::setNeutralLossMonoMass(DoubleReal mono_mass)
  {
    neutral_loss_mono_mass_ = mono_mass;
  }

  DoubleReal ResidueModification::getNeutralLossMonoMass() const
  {
    return neutral_loss_mono_mass_;
  }

  void ResidueModification::setNeutralLossAverageMass(DoubleReal average_mass)
  {
    neutral_loss_average_mass_ = average_mass;
  }

  DoubleReal ResidueModification::getNeutralLossAverageMass() const
  {
    return neutral_loss_average_mass_;
  }

  bool ResidueModification::hasNeutralLoss() const
  {
    return neutral_loss_diff_formula_ != "";
  }

}
