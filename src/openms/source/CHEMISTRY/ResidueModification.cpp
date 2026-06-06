// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>
#include <iostream>
#include <utility>

using namespace std;

namespace OpenMS
{

  ResidueModification::ResidueModification() :
    unimod_record_id_(-1),
    term_spec_(ResidueModification::ANYWHERE),
    origin_('X'),
    classification_(ResidueModification::ARTIFACT),
    average_mass_(0.0),
    mono_mass_(0.0),
    diff_average_mass_(0.0),
    diff_mono_mass_(0.0),
    neutral_loss_mono_masses_(0),
    neutral_loss_average_masses_(0)
  {
  }

  bool ResidueModification::operator<(const ResidueModification& rhs) const
  {
    return std::tie(
           id_,
           full_id_,
           psi_mod_accession_,
           unimod_record_id_,
           full_name_,
           name_,
           term_spec_,
           origin_,
           classification_,
           average_mass_,
           mono_mass_,
           diff_average_mass_,
           diff_mono_mass_,
           formula_,
           diff_formula_,
           synonyms_,
           neutral_loss_diff_formulas_,
           neutral_loss_mono_masses_,
           neutral_loss_average_masses_
    ) < std::tie(
           rhs.id_,
           rhs.full_id_,
           rhs.psi_mod_accession_,
           rhs.unimod_record_id_,
           rhs.full_name_,
           rhs.name_,
           rhs.term_spec_,
           rhs.origin_,
           rhs.classification_,
           rhs.average_mass_,
           rhs.mono_mass_,
           rhs.diff_average_mass_,
           rhs.diff_mono_mass_,
           rhs.formula_,
           rhs.diff_formula_,
           rhs.synonyms_,
           rhs.neutral_loss_diff_formulas_,
           rhs.neutral_loss_mono_masses_,
           rhs.neutral_loss_average_masses_
    );
  }


  bool ResidueModification::operator==(const ResidueModification& rhs) const
  {
    return id_ == rhs.id_ &&
           full_id_ == rhs.full_id_ &&
           psi_mod_accession_ == rhs.psi_mod_accession_ &&
           unimod_record_id_ == rhs.unimod_record_id_ &&
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
           neutral_loss_diff_formulas_ == rhs.neutral_loss_diff_formulas_ &&
           neutral_loss_mono_masses_ == rhs.neutral_loss_mono_masses_ &&
           neutral_loss_average_masses_ == rhs.neutral_loss_average_masses_;
  }

  bool ResidueModification::operator!=(const ResidueModification& rhs) const
  {
    return !(*this == rhs);
  }

  ResidueModification::~ResidueModification() = default;

  void ResidueModification::setId(const std::string& id)
  {
    id_ = id;
  }

  const std::string& ResidueModification::getId() const
  {
    return id_;
  }

  void ResidueModification::setFullId(const std::string& full_id)
  {
    if (full_id.empty())
    {
      if (id_.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot create full ID for modification with missing (short) ID.");
      }
      std::string specificity;
      if (term_spec_ != ResidueModification::ANYWHERE)
      {
        specificity = getTermSpecificityName(); // "C-term" or "N-term"
      }
      if (!specificity.empty() && (origin_ != 'X'))
      {
        specificity += " " + StringUtils::toStr(origin_);
      }
      else if (specificity.empty())
      {
        specificity = origin_; // shouldn't be "X" in this case
      }
      full_id_ = id_ + " (" + specificity + ")";
    }
    else
    {
      full_id_ = full_id;
    }
  }

  const std::string& ResidueModification::getFullId() const
  {
    return full_id_;
  }

  void ResidueModification::setPSIMODAccession(const std::string& id)
  {
    psi_mod_accession_ = id;
  }

  const std::string& ResidueModification::getPSIMODAccession() const
  {
    return psi_mod_accession_;
  }

  void ResidueModification::setUniModRecordId(const Int& id)
  {
    unimod_record_id_ = id;
  }

  const Int& ResidueModification::getUniModRecordId() const
  {
    return unimod_record_id_;
  }

  const std::string ResidueModification::getUniModAccession() const
  {
    if (unimod_record_id_ < 0) return "";
    return StringUtils::toStr("UniMod:") + unimod_record_id_; // return copy of temp object
  }

  void ResidueModification::setFullName(const std::string& full_name)
  {
    full_name_ = full_name;
  }

  const std::string& ResidueModification::getFullName() const
  {
    return full_name_;
  }

  void ResidueModification::setName(const std::string& name)
  {
    name_ = name;
  }

  const std::string& ResidueModification::getName() const
  {
    return name_;
  }

  void ResidueModification::setTermSpecificity(TermSpecificity term_spec)
  {
    if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not a valid terminal specificity",StringUtils::toStr(term_spec));
    }
    term_spec_ = term_spec;
  }

  void ResidueModification::setTermSpecificity(const std::string& term_spec)
  {
    if (term_spec == "C-term")
    {
      term_spec_ = C_TERM;
    }
    else if (term_spec == "N-term")
    {
      term_spec_ = N_TERM;
    }
    else if (term_spec == "none")
    {
      term_spec_ = ANYWHERE;
    }
    else if (term_spec == "Protein N-term")
    {
      term_spec_ = PROTEIN_N_TERM;
    }
    else if (term_spec == "Protein C-term")
    {
      term_spec_ = PROTEIN_C_TERM;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not a valid terminal specificity", term_spec);
    }
  }

  ResidueModification::TermSpecificity ResidueModification::getTermSpecificity() const
  {
    return term_spec_;
  }

  std::string ResidueModification::getTermSpecificityName(TermSpecificity term_spec) const
  {
    if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
    {
      term_spec = term_spec_;
    }
    switch (term_spec)
    {
    case C_TERM: return "C-term";

    case N_TERM: return "N-term";

    case PROTEIN_C_TERM: return "Protein C-term";

    case PROTEIN_N_TERM: return "Protein N-term";

    case ANYWHERE: return "none";

    default: break; // shouldn't happen
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No name for this terminal specificity",StringUtils::toStr(term_spec));
  }

  void ResidueModification::setOrigin(char origin)
  {
    if ((origin >= 'A') && (origin <= 'Y') && (origin != 'B') && (origin != 'J'))
    {
      origin_ = origin;
    }
    else if ((origin >= 'a') && (origin <= 'y') && (origin != 'b') && (origin != 'j'))
    {
      origin_ = toupper(origin);
    }
    else
    {
      std::string msg = "Modification '" + id_ + "': origin must be a letter from A to Y, excluding B and J.";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg,StringUtils::toStr(origin));
    }
  }

  char ResidueModification::getOrigin() const
  {
    return origin_;
  }

  void ResidueModification::setSourceClassification(SourceClassification classification)
  {
    classification_ = classification;
  }

  void ResidueModification::setSourceClassification(const std::string& classification)
  {
    std::string c = classification;
    StringUtils::toLower(c);
    if (c == "artifact" || c == "artefact") // unimod uses Artefact (BE) not Artifact (AE)
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

  ResidueModification::SourceClassification ResidueModification::getSourceClassification() const
  {
    return classification_;
  }

  std::string ResidueModification::getSourceClassificationName(SourceClassification classification) const
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
  }

  void ResidueModification::setAverageMass(double mass)
  {
    average_mass_ = mass;
  }

  double ResidueModification::getAverageMass() const
  {
    return average_mass_;
  }

  void ResidueModification::setMonoMass(double mass)
  {
    mono_mass_ = mass;
  }

  double ResidueModification::getMonoMass() const
  {
    return mono_mass_;
  }

  void ResidueModification::setDiffAverageMass(double mass)
  {
    diff_average_mass_ = mass;
  }

  double ResidueModification::getDiffAverageMass() const
  {
    return diff_average_mass_;
  }

  void ResidueModification::setDiffMonoMass(double mass)
  {
    diff_mono_mass_ = mass;
  }

  double ResidueModification::getDiffMonoMass() const
  {
    return diff_mono_mass_;
  }

  void ResidueModification::setFormula(const std::string& formula)
  {
    formula_ = formula;
  }

  const std::string& ResidueModification::getFormula() const
  {
    return formula_;
  }

  void ResidueModification::setDiffFormula(const EmpiricalFormula& diff_formula)
  {
    diff_formula_ = diff_formula;
  }

  const EmpiricalFormula& ResidueModification::getDiffFormula() const
  {
    return diff_formula_;
  }

  void ResidueModification::addSynonym(const std::string& synonym)
  {
    synonyms_.insert(synonym);
  }

  void ResidueModification::setSynonyms(const set<std::string>& synonyms)
  {
    synonyms_ = synonyms;
  }

  const set<std::string>& ResidueModification::getSynonyms() const
  {
    return synonyms_;
  }

  void ResidueModification::setNeutralLossDiffFormulas(const vector<EmpiricalFormula>& diff_formulas)
  {
    neutral_loss_diff_formulas_ = diff_formulas;
  }

  const vector<EmpiricalFormula>& ResidueModification::getNeutralLossDiffFormulas() const
  {
    return neutral_loss_diff_formulas_;
  }

  void ResidueModification::setNeutralLossMonoMasses(vector<double> mono_masses)
  {
    neutral_loss_mono_masses_ = std::move(mono_masses);
  }

  vector<double> ResidueModification::getNeutralLossMonoMasses() const
  {
    return neutral_loss_mono_masses_;
  }

  void ResidueModification::setNeutralLossAverageMasses(vector<double> average_masses)
  {
    neutral_loss_average_masses_ = std::move(average_masses);
  }

  vector<double> ResidueModification::getNeutralLossAverageMasses() const
  {
    return neutral_loss_average_masses_;
  }

  bool ResidueModification::hasNeutralLoss() const
  {
    return !neutral_loss_diff_formulas_.empty() && !neutral_loss_diff_formulas_[0].isCharged();
  }

  bool ResidueModification::isUserDefined() const
  {
    return id_.empty() && !full_id_.empty();
  }
  
  const ResidueModification* ResidueModification::createUnknownFromMassString(const std::string& mod,
                                                                              const double mass,
                                                                              const bool delta_mass,
                                                                              const TermSpecificity specificity,
                                                                              const Residue* residue)
  {
    const ModificationsDB* mod_db = ModificationsDB::getInstance();

    // -----------------------------------
    // Dealing with an unknown modification
    // -----------------------------------

    // Notes on mass calculation: AASequence::getMonoWeight uses DiffMonoMass
    // for its calculation of C/N-terminal modification mass and it uses
    // getMonoWeight(Residue::Internal) for each Residue. The Residue weight is
    // set when adding a modification using setModification_
    if (specificity == ResidueModification::N_TERM || specificity == ResidueModification::PROTEIN_N_TERM)
    {
      std::string residue_name = "[" + mod + "]";
      std::string residue_id = ".n" + residue_name;

      // Check if it already exists, if not create new modification, transfer
      // ownership to ModDB
      if (!mod_db->has(residue_id))
      {
        unique_ptr<ResidueModification> new_mod(new ResidueModification);
        new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
        new_mod->setFullName(residue_name); // display name
        new_mod->setTermSpecificity(specificity);

        // set masses
        if (delta_mass)
        {
          new_mod->setMonoMass(mass + Residue::getInternalToNTerm().getMonoWeight());
          // new_mod->setAverageMass(mass + residue->getAverageWeight());
          new_mod->setDiffMonoMass(mass);
        }
        else
        {
          new_mod->setMonoMass(mass);
          // new_mod->setAverageMass(mass);
          new_mod->setDiffMonoMass(mass - Residue::getInternalToNTerm().getMonoWeight());
        }
     
        return mod_db->addModification(std::move(new_mod));
      }
      else
      {
        Size mod_idx = mod_db->findModificationIndex(residue_id);
        return mod_db->getModification(mod_idx);
      }
    }
    else
      if (specificity == ResidueModification::C_TERM || specificity == ResidueModification::PROTEIN_C_TERM)
      {
        std::string residue_name = "[" + mod + "]";
        std::string residue_id = ".c" + residue_name;

        // Check if it already exists, if not create new modification, transfer
        // ownership to ModDB
        if (!mod_db->has(residue_id))
        {
          unique_ptr<ResidueModification> new_mod(new ResidueModification);
          new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
          new_mod->setFullName(residue_name); // display name
          new_mod->setTermSpecificity(specificity);

          // set masses
          if (delta_mass)
          {
            new_mod->setMonoMass(mass + Residue::getInternalToCTerm().getMonoWeight());
            // new_mod->setAverageMass(mass + residue->getAverageWeight());
            new_mod->setDiffMonoMass(mass);
          }
          else
          {
            new_mod->setMonoMass(mass);
            // new_mod->setAverageMass(mass);
            new_mod->setDiffMonoMass(mass - Residue::getInternalToCTerm().getMonoWeight());
          }

          return mod_db->addModification(std::move(new_mod));
        }
        else
        {
          Size mod_idx = mod_db->findModificationIndex(residue_id);
          return mod_db->getModification(mod_idx);
        }
      }
      else
      {
        if (residue == nullptr)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot create non-terminal mod without origin AA residue.", "nullptr");
        }
        std::string modification_name = "[" + mod + "]";
        std::string residue_id =StringUtils::toStr(residue->getOneLetterCode()) + modification_name; // e.g. N[12345.6]

        if (!mod_db->has(residue_id))
        {
          // create new modification
          unique_ptr<ResidueModification> new_mod(new ResidueModification);
          new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
          new_mod->setFullName(modification_name); // display name

          // We will set origin to make sure the same modification will be used
          // for the same AA
          new_mod->setOrigin(residue->getOneLetterCode()[0]);

          // set masses
          if (delta_mass)
          {
            new_mod->setMonoMass(mass + residue->getMonoWeight());
            new_mod->setAverageMass(mass + residue->getAverageWeight());
            new_mod->setDiffMonoMass(mass);
          }
          else
          {
            new_mod->setMonoMass(mass);
            new_mod->setAverageMass(mass);
            new_mod->setDiffMonoMass(mass - residue->getMonoWeight());
          }

          return mod_db->addModification(std::move(new_mod));
        }
        else
        {
          Size mod_idx = mod_db->findModificationIndex(residue_id);
          return mod_db->getModification(mod_idx);
        }
      }
  }
  
  const ResidueModification* ResidueModification::combineMods(const ResidueModification* base, 
                                                              const std::set<const ResidueModification*>& addons, 
                                                              bool allow_unknown_masses,
                                                              const Residue* residue)
  {
    if (addons.empty() && base == nullptr)
    {
      return nullptr;
    }
    if (base != nullptr && base->isUserDefined() && !allow_unknown_masses)
    { // do not apply anything to an already merged mod
      OPENMS_LOG_INFO << "Note: Invalid merge operation on already merged/user-defined modification!\n";
      return base;
    }
    
    const ResidueModification* mod_merged = base;

    std::set<const ResidueModification*>::const_iterator it_addon = addons.begin();

    if (base == nullptr) // base is empty
    {
      mod_merged = *it_addon;
      ++it_addon; // start adding at the second element
    }
    if (it_addon == addons.end()) return mod_merged; // only one mod in total given
    
    // compute combined mass ...    
    double new_mass{ mod_merged->getDiffMonoMass() };
    for (const auto& mod_new : addons)
    { // build a new mod, combining the existing and the new one
      if (mod_merged->getTermSpecificity() != mod_new->getTermSpecificity())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modifications to be merged to not have the same term specificity: " + mod_merged->getTermSpecificityName() + " != " + mod_new->getTermSpecificityName());
      }
      if (mod_merged->getOrigin() != mod_new->getOrigin())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Modifications to be merged to not have the same origin: ") + mod_merged->getOrigin() + " != " + mod_new->getOrigin());
      }
      new_mass += mod_new->getDiffMonoMass();
    }

    // sanity check: mods and residue need same origin
    if (mod_merged->getTermSpecificity() == ANYWHERE && residue != nullptr && residue->getOneLetterCode()[0] != mod_merged->getOrigin())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Modification and residue do not have the same origin: ") + mod_merged->getOrigin() + " != " + residue->getOneLetterCode());
    }
    // create new residue from it
    const ResidueModification* mod_sum =
      ResidueModification::createUnknownFromMassString(getDiffMonoMassString(new_mass),
                                                       new_mass,
                                                       true,
                                                       mod_merged->getTermSpecificity(),
                                                       residue);

    return mod_sum;
  }

  std::string ResidueModification::toString() const
  {
    std::string ret;

    if (term_spec_ != ANYWHERE) ret = ".";
    else ret = origin_;

    if (!id_.empty())
    {
      ret.reserve(id_.size() + 3);
      ret += "(";
      ret += id_;
      ret += ")";
      return ret;
    }
    
    if (!getFullName().empty())
    { // using this is questionable ... it will not work for default mods (which do have id_ defined, so this case will never happen... but still)
      ret += getFullName();
      return ret;
    }
    
    if (isUserDefined()) // id_ empty but full_id_ is present
    { // user-defined modification
      if (diff_mono_mass_ != 0.0)
      {
        ret += getDiffMonoMassWithBracket(diff_mono_mass_);
      }
      else if (mono_mass_ != 0.0)
      {
        ret += getMonoMassWithBracket(mono_mass_);
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Residue has an invalid user-defined modification. This is a bug. Please report it!", "");
      }
      return ret;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ResidueModification is in an invalid state. This is a bug. Please report it!", "");
  }
  std::string ResidueModification::getDiffMonoMassString(const double diff_mono_mass)
  {
    return StringUtils::toStr(diff_mono_mass < 0.0 ? "-" : "+") += std::fabs(diff_mono_mass);
  }
  std::string ResidueModification::getDiffMonoMassWithBracket(const double diff_mono_mass)
  {
    std::string ret(1, '[');
    ret += getDiffMonoMassString(diff_mono_mass);
    ret += ']';
    return ret;
  }
  std::string ResidueModification::getMonoMassWithBracket(const double mono_mass)
  {
    if (mono_mass < 0.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Modification has negative mono mass. Cannot distinguish between delta masses due to '-'!",StringUtils::toStr(mono_mass));
    }
    std::string ret(1, '[');
    ret += mono_mass;
    ret += ']';
    return ret;
  }
  
}
