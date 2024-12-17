// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <cmath>
#include <algorithm>
#include <map>

using namespace std;

namespace OpenMS
{

  const ResidueModification* proteinTerminalResidueHelper( ModificationsDB* mod_db,
      const char term,
      const std::string& str,
      const std::string& mod,
      const std::string& res)
  {
    ResidueModification::TermSpecificity protein_term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;

    if (term == 'c')
    {
      protein_term_spec = ResidueModification::PROTEIN_C_TERM;
    }
    else if (term == 'n')
    {
      protein_term_spec = ResidueModification::PROTEIN_N_TERM;
    }

    const ResidueModification* m;
    try
    {
      m = mod_db->getModification(mod, res, protein_term_spec);
    }
    catch (...)
    {
      // catch and rethrow with some additional information on term
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
          "Cannot convert string to peptide modification. No terminal modification matches to term specificity and origin.");
    }
    return m;
  }

  const ResidueModification* terminalResidueHelper( ModificationsDB* mod_db,
      const char term,
      bool protein_term,
      const std::string& str,
      const std::string& mod,
      const std::string& res)
  {
    ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;
    ResidueModification::TermSpecificity protein_term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;;

    if (term == 'c')
    {
      term_spec = ResidueModification::C_TERM;
      protein_term_spec = ResidueModification::PROTEIN_C_TERM;
    }
    else if (term == 'n')
    {
      term_spec = ResidueModification::N_TERM;
      protein_term_spec = ResidueModification::PROTEIN_N_TERM;
    }

    const ResidueModification* m;
    try
    {
      m = mod_db->getModification(mod, res, term_spec);
    }
    catch (...)
    {
      if (!protein_term)
      {
        // catch and rethrow with some additional information on term
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
            "Cannot convert string to peptide modification. No terminal modification matches to term specificity and origin.");
      }

      try
      {
        m = mod_db->getModification(mod, res, protein_term_spec);
      }
      catch (...)
      {
        // catch and rethrow with some additional information on term
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
            "Cannot convert string to peptide modification. No terminal modification matches to term specificity and origin.");
      }
    }
    return m;
  }

  AASequence::AASequence(const String& s)
  {
    parseString_(s, *this, true);
  }

  AASequence::AASequence(const char* s)
  {
    parseString_(s, *this, true);
  }

  AASequence::AASequence(const String& s, bool permissive)
  {
    parseString_(s, *this, permissive);
  }

  AASequence::AASequence(const char* s, bool permissive)
  {
    parseString_(s, *this, permissive);
  }

  const Residue& AASequence::getResidue(Size index) const
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, peptide_.size());
    }
    return *peptide_[index];
  }

  String AASequence::toString() const
  {
    std::stringstream ss;
    ss << *this;
    return String(ss.str());
  }

  String AASequence::toUnmodifiedString() const
  {
    String tmp;
    for (ConstIterator it = begin(); it != end(); ++it)
    {
      tmp += it->getOneLetterCode();
    }
    return tmp;
  }

  String AASequence::toUniModString() const
  {
    const AASequence & seq = *this;

    String bs;
    if (seq.empty()) return bs;

    if (seq.hasNTerminalModification())
    {
      const ResidueModification& mod = *(seq.getNTerminalModification());
      double nominal_mass = Residue::getInternalToNTerm().getMonoWeight() + mod.getDiffMonoMass();

      if (mod.getUniModRecordId() > -1)
      {
        bs += ".(" + mod.getUniModAccession() + ")";
      }
      else
      {
        bs += ".[" + String(nominal_mass) + "]";
      }
    }

    for (Size i = 0; i != seq.size(); ++i)
    {
      const Residue& r = seq[i];
      const String& aa = r.getOneLetterCode();
      if (r.isModified())
      {
        const ResidueModification& mod = *(r.getModification());
        double nominal_mass = r.getMonoWeight(Residue::Internal);

        if (mod.getUniModRecordId() > -1)
        {
          bs += aa + "(" + mod.getUniModAccession() + ")";
        }
        else
        {
          bs += aa + "[" + String(nominal_mass) + "]";
        }
      }
      else  // amino acid not modified
      {
        bs += aa;
      }
    }

    if (seq.hasCTerminalModification())
    {
      const ResidueModification& mod = *(seq.getCTerminalModification());
      double nominal_mass = Residue::getInternalToCTerm().getMonoWeight() + mod.getDiffMonoMass();

      if (mod.getUniModRecordId() > -1)
      {
        bs += ".(" + mod.getUniModAccession() + ")";
      }
      else
      {
        bs += ".[" + String(nominal_mass) + "]";
      }
    }
    return bs;
  }

  String AASequence::toBracketString(bool integer_mass, bool mass_delta, const vector<String> & fixed_modifications) const
  {
    const AASequence& seq = *this;

    String bs;
    if (seq.empty())
    {
      return bs;
    }
    if (seq.hasNTerminalModification())
    {
      const ResidueModification& mod = *(seq.getNTerminalModification());
      const String & nterm_mod_name = mod.getFullId(); // e.g. "Acetyl (N-term)"

      // only add to string if not a fixed modification
      if (std::find(fixed_modifications.begin(), fixed_modifications.end(), nterm_mod_name) == fixed_modifications.end())
      {
        double nominal_mass = mod.getDiffMonoMass();
       if (!mass_delta)
        {
          nominal_mass += Residue::getInternalToNTerm().getMonoWeight();
        }
        String sign = (mass_delta && nominal_mass > 0) ? "+" : ""; // the '-' will be printed automatically by conversion to string
        if (integer_mass)
        {
          bs += String("n[") + sign + static_cast<int>(std::round(nominal_mass)) + "]";
        }
        else
        {
          bs += "n[" + sign + String(nominal_mass) + "]";
        }
      }
    }

    for (Size i = 0; i != seq.size(); ++i)
    {
      const Residue& r = seq[i];
      const String aa = !r.getOneLetterCode().empty() ? r.getOneLetterCode() : "X";
      if (r.isModified())
      {
        const ResidueModification& mod = *(r.getModification());

        const String & mod_name = mod.getFullId();
        if (std::find(fixed_modifications.begin(), fixed_modifications.end(), mod_name) == fixed_modifications.end())
        {
          double nominal_mass;
          if (mass_delta)
          {
            nominal_mass = mod.getDiffMonoMass();
          }
          else
          {
            nominal_mass = r.getMonoWeight(Residue::Internal);
          }
          String sign = (mass_delta && nominal_mass > 0) ? "+" : "";
          if (aa == "X")
          {
            // cannot have delta mass for X
            nominal_mass = r.getMonoWeight(Residue::Internal);
            sign = "";
          }

          if (integer_mass)
          {
            bs += aa + String("[") + sign + static_cast<int>(std::round(nominal_mass)) + "]";
          }
          else
          {
            bs += aa + "[" + sign + String(nominal_mass) + "]";
          }
        }
        else
        {
          bs += aa; // don't print fixed modification
        }
      }
      else  // amino acid not modified
      {
        bs += aa;
      }
    }

    if (seq.hasCTerminalModification())
    {
      const ResidueModification& mod = *(seq.getCTerminalModification());
      const String & cterm_mod_name = mod.getFullId();

      // only add to string if not a fixed modification
      if (std::find(fixed_modifications.begin(), fixed_modifications.end(), cterm_mod_name) == fixed_modifications.end())
      {
        double nominal_mass = mod.getDiffMonoMass();
        if (!mass_delta)
        {
          nominal_mass += Residue::getInternalToCTerm().getMonoWeight();
        }
        String sign = (mass_delta && nominal_mass > 0) ? "+" : "";
        if (integer_mass)
        {
          bs += String("c[") + sign + static_cast<int>(std::round(nominal_mass)) + "]";
        }
        else
        {
          bs += "c[" + sign + String(nominal_mass) + "]";
        }
      }
    }
    return bs;
  }

  bool AASequence::operator<(const AASequence& rhs) const
  {
    // check size
    if (peptide_.size() != rhs.peptide_.size())
    {
      return (peptide_.size() < rhs.peptide_.size());
    }

    // when checking terminal mods, "no mod" is less than "any mod"
    if (n_term_mod_ && !rhs.n_term_mod_)
    {
      return false;
    }
    else if (!n_term_mod_ && rhs.n_term_mod_)
    {
      return true;
    }
    else if (n_term_mod_ && rhs.n_term_mod_ && (n_term_mod_ != rhs.n_term_mod_))
    {
      return (n_term_mod_->getId() < rhs.n_term_mod_->getId());
    }

    ConstIterator a = begin();
    ConstIterator b = rhs.begin();

    // check one letter codes
    for (; a != end(); ++a, ++b)
    {
      if (a->getOneLetterCode() != b->getOneLetterCode())
      {
        return (a->getOneLetterCode() < b->getOneLetterCode());
      }
      else if (a->getModification() != b->getModification())
      {
        return (a->getModification() < b->getModification());
      }
    }

    // c-term
    if (c_term_mod_ && !rhs.c_term_mod_)
    {
      return false;
    }
    else if (!c_term_mod_ && rhs.c_term_mod_)
    {
      return true;
    }
    else if (c_term_mod_ && rhs.c_term_mod_ && (c_term_mod_ != rhs.c_term_mod_))
    {
      return (c_term_mod_->getId() < rhs.c_term_mod_->getId());
    }

    return false;
  }

  EmpiricalFormula AASequence::getFormula(Residue::ResidueType type, Int charge) const
  {
    if (!peptide_.empty())
    {
      // Initialize with the missing/additional protons
      EmpiricalFormula ef; // = EmpiricalFormula("H") * charge; ??
      ef.setCharge(charge);

      // terminal modifications
      if (n_term_mod_ != nullptr &&
        (type == Residue::Full || type == Residue::AIon ||
         type == Residue::BIon || type == Residue::CIon ||
         type == Residue::NTerminal))
      {
        ef += n_term_mod_->getDiffFormula();
      }


      if (c_term_mod_ != nullptr &&
        (type == Residue::Full || type == Residue::XIon ||
         type == Residue::YIon || type == Residue::ZIon ||
         type == Residue::CTerminal))
      {
        ef += c_term_mod_->getDiffFormula();
      }

      static auto const rx = ResidueDB::getInstance()->getResidue("X");
      for (auto const& e : peptide_)
      {
        // While PEPTIX[123]DE makes sense and represents an unknown mass of 123.0
        // Da, the sequence PEPTIXDE does not make sense as it is unclear what a
        // standard internal residue including named modifications
        if (e == rx)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot get EF of sequence with unknown AA 'X'.", toString());
        }
        ef += e->getFormula(Residue::Internal);
      }

      // add the missing formula part
      switch (type)
      {
        case Residue::Full:
        {
          return ef + Residue::getInternalToFull();
        }
        case Residue::Internal:
        {
          return ef;
        }
        case Residue::NTerminal:
        {
          return ef + Residue::getInternalToNTerm();
        }
        case Residue::CTerminal:
        {
          return ef + Residue::getInternalToCTerm();
        }
        case Residue::AIon:
        {
          return ef + Residue::getInternalToAIon();
        }
        case Residue::BIon:
        {
          return ef + Residue::getInternalToBIon();
        }
        case Residue::CIon:
        {
          return ef + Residue::getInternalToCIon();
        }
        case Residue::XIon:
        {
          return ef + Residue::getInternalToXIon();
        }
        case Residue::YIon:
        {
          return ef + Residue::getInternalToYIon();
        }
        case Residue::ZIon:
        {
          return ef + Residue::getInternalToZIon();
        }
        default:
          OPENMS_LOG_ERROR << "AASequence::getFormula: unknown ResidueType" << std::endl;
      }

      return ef;
    }
    else
    {
      OPENMS_LOG_ERROR << "AASequence::getFormula: Formula for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
      return EmpiricalFormula("");
    }
  }

  double AASequence::getAverageWeight(Residue::ResidueType type, Int charge) const
  {
    // check whether tags are present
    double tag_offset(0);
    static auto const rx = ResidueDB::getInstance()->getResidue("X");
    for (auto const& e : peptide_)
    {
      // While PEPTIX[123]DE makes sense and represents an unknown mass of 123.0
      // Da, the sequence PEPTIXDE does not make sense as it is unclear what a
      // standard internal residue including named modifications
      if (e == rx) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot get weight of sequence with unknown AA 'X' with unknown mass.", toString());
      if (e->getOneLetterCode().empty())
      {
        tag_offset += e->getAverageWeight(Residue::Internal);
      }
    }
    // TODO inefficient, if averageWeight is already set in the Residue
    return tag_offset + getFormula(type, charge).getAverageWeight();
  }

  double AASequence::getMZ(Int charge, Residue::ResidueType type) const
  {
    if (charge == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can't calculate mass-to-charge ratio for charge=0.", toString());
    }
    return getMonoWeight(type, charge) / charge;
  }

  double AASequence::getMonoWeight(Residue::ResidueType type, Int charge) const
  {
    if (!peptide_.empty())
    {
      double mono_weight(Constants::PROTON_MASS_U * charge);

      // terminal modifications
      if (n_term_mod_ != nullptr &&
          (type == Residue::Full || type == Residue::AIon ||
           type == Residue::BIon || type == Residue::CIon ||
           type == Residue::NTerminal))
      {
        mono_weight += n_term_mod_->getDiffMonoMass();
      }

      if (c_term_mod_ != nullptr &&
          (type == Residue::Full || type == Residue::XIon ||
           type == Residue::YIon || type == Residue::ZIon ||
           type == Residue::CTerminal))
      {
        mono_weight += c_term_mod_->getDiffMonoMass();
      }
      static auto const rx = ResidueDB::getInstance()->getResidue("X");
      for (auto const& e : peptide_)
      {
        // While PEPTIX[123]DE makes sense and represents an unknown mass of 123.0
        // Da, the sequence PEPTIXDE does not make sense as it is unclear what a
        // standard internal residue including named modifications
        if (e == rx)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot get weight of sequence with unknown AA 'X' with unknown mass.", toString());
        }
        // single, unknown residue should represent.
        mono_weight += e->getMonoWeight(Residue::Internal);
      }

      // add the missing formula part
      switch (type)
      {
        case Residue::Full:
        {
          return mono_weight + Residue::getInternalToFull().getMonoWeight();
        }
        case Residue::Internal:
        {
          return mono_weight;
        }
        case Residue::NTerminal:
        {
          return mono_weight + Residue::getInternalToNTerm().getMonoWeight();
        }
        case Residue::CTerminal:
        {
          return mono_weight + Residue::getInternalToCTerm().getMonoWeight();
        }
        case Residue::AIon:
        {
          return mono_weight + Residue::getInternalToAIon().getMonoWeight();
        }
        case Residue::BIon:
        {
          return mono_weight + Residue::getInternalToBIon().getMonoWeight();
        }
        case Residue::CIon:
        {
          return mono_weight + Residue::getInternalToCIon().getMonoWeight();
        }
        case Residue::XIon:
        {
          return mono_weight + Residue::getInternalToXIon().getMonoWeight();
        }
        case Residue::YIon:
        {
          return mono_weight + Residue::getInternalToYIon().getMonoWeight();
        }
        case Residue::ZIon:
        {
          return mono_weight + Residue::getInternalToZIon().getMonoWeight();
        }
        default:
          OPENMS_LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
      }

      return mono_weight;
    }
    else
    {
      OPENMS_LOG_ERROR << "AASequence::getMonoWeight: Mass for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
      return 0.0;
    }
}




/*void AASequence::getNeutralLosses(Map<const EmpiricalFormula, UInt) const
  {
      // the following losses are from the Zhang paper (AC, 76, 14, 2004)
      // charge directed*/
/*
  static const EmpiricalFormula R_44("NH2CHNH");
  static const EmpiricalFormula R_59("CN3H5"); // guanidine
  static const EmpiricalFormula R_61("N2H4CH");
  // charge remote
  static const EmpiricalFormula R_60("N2H4CO"); // combination of NH=C=NH + C-terminal H2O
  static const EmpiricalFormula H2O("H2O"); // loss from the C-terminus
  static const EmpiricalFormula NH3("NH3");
  Map<const EmpiricalFormula*, UInt> losses;

  for (Size i=0;i!=peptide_.size();++i)
  {
      if (peptide_[i]->hasNeutralLoss())
      {
          const EmpiricalFormula* loss = peptide_[i]->getLossFormulas();
          if (losses.find(loss) != losses.end())
          {
              losses[loss]++;
          }
          else
          {
              losses[loss] = 1;
          }
      }


      // TODO: hack this should be in the data file
      if (peptide_[i]->getOneLetterCode() == "R")
      {
          losses[&R_44] = 1;
          losses[&R_59] = 1;
          losses[&R_61] = 1;
          losses[&R_60] = 1;
      }
      losses[&H2O] = 1;
      losses[&NH3] = 1;
  }
  return losses;
}*/

  const Residue& AASequence::operator[](Size index) const
  {
    if (index >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, size());
    }
    return *peptide_[index];
  }

  AASequence& AASequence::operator+=(const AASequence& sequence)
  {
    for (Size i = 0; i != sequence.peptide_.size(); ++i)
    {
      peptide_.push_back(sequence.peptide_[i]);
    }
    return *this;
  }

  AASequence AASequence::operator+(const AASequence& sequence) const
  {
    AASequence seq;
    seq.peptide_ = peptide_;
    for (Size i = 0; i != sequence.peptide_.size(); ++i)
    {
      seq.peptide_.push_back(sequence.peptide_[i]);
    }
    return seq;
  }

  AASequence AASequence::operator+(const Residue* residue) const
  {
    if (!ResidueDB::getInstance()->hasResidue(residue))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "given residue");
    }
    AASequence seq = *this;
    seq += residue;
    return seq;
  }

  AASequence& AASequence::operator+=(const Residue* residue)
  {
    if (!ResidueDB::getInstance()->hasResidue(residue))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "given residue");
    }
    peptide_.push_back(residue);
    return *this;
  }

  Size AASequence::size() const
  {
    return peptide_.size();
  }

  AASequence AASequence::getPrefix(Size index) const
  {
    if (index > size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, size());
    }
    if (index == size())
    {
      return *this;
    }

    AASequence seq;
    seq.n_term_mod_ = n_term_mod_;
    seq.peptide_.reserve(index);
    seq.peptide_.insert(seq.peptide_.end(), peptide_.begin(), peptide_.begin() + index);
    return seq;
  }

  AASequence AASequence::getSuffix(Size index) const
  {
    if (index > size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, size());
    }

    if (index == size())
    {
      return *this;
    }

    AASequence seq;
    seq.c_term_mod_ = c_term_mod_;
    seq.peptide_.reserve(size() - index);
    seq.peptide_.insert(seq.peptide_.end(), peptide_.begin() + (size() - index), peptide_.end());
    return seq;
  }

  AASequence AASequence::getSubsequence(Size index, UInt num) const
  {
    if (index >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, size());
    }
    if (index + num > size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index + num, size());
    }

    AASequence seq;
    if (index == 0)
    {
      seq.n_term_mod_ = n_term_mod_;
    }
    if (index + num == this->size())
    {
      seq.c_term_mod_ = c_term_mod_;
    }
    seq.peptide_.reserve(num);
    seq.peptide_.insert(seq.peptide_.end(), peptide_.begin() + index, peptide_.begin() + index + num);

    return seq;
  }

  bool AASequence::has(const Residue& residue) const
  {
    for (const Residue* rp : peptide_)
    {
      if (*rp == residue) return true;
    }
    return false;
  }

  bool AASequence::hasSubsequence(const AASequence& sub) const
  {
    if (sub.empty())
    {
      return true;
    }
    Size full_size = peptide_.size();
    Size sub_size = sub.peptide_.size();
    if (sub_size > full_size)
    {
      return false;
    }

    for (Size i = 0; i < full_size - sub_size + 1; ++i)
    {
      if (peptide_[i] == sub.peptide_[0])
      {
        Size j = 1;
        for (; j < sub_size; ++j)
        {
          if (peptide_[j + i] != sub.peptide_[j])
          {
            break;
          }
        }
        // check if we reached last position
        if (j == sub_size)
        {
          return true;
        }
      }
    }
    return false;
  }

  bool AASequence::hasPrefix(const AASequence& sequence) const
  {
    if (sequence.empty())
    {
      return true;
    }
    if (sequence.size() > peptide_.size())
    {
      return false;
    }
    if (sequence.n_term_mod_ != n_term_mod_)
    {
      return false;
    }
    if (sequence.size() == peptide_.size() && sequence.c_term_mod_ != c_term_mod_)
    {
      return false;
    }

    for (Size i = 0; i != sequence.size(); ++i)
    {
      if (sequence.peptide_[i] != peptide_[i])
      {
        return false;
      }
    }
    return true;
  }

  bool AASequence::hasSuffix(const AASequence& sequence) const
  {
    if (sequence.empty())
    {
      return true;
    }
    if (sequence.size() > peptide_.size())
    {
      return false;
    }
    if (sequence.c_term_mod_ != c_term_mod_)
    {
      return false;
    }
    if (sequence.size() == peptide_.size() && sequence.n_term_mod_ != n_term_mod_)
    {
      return false;
    }

    for (Size i = 0; i != sequence.size(); ++i)
    {
      if (sequence.peptide_[sequence.size() - 1 - i] != peptide_[size() - 1 - i])
      {
        return false;
      }
    }
    return true;
  }

  bool AASequence::operator==(const AASequence& peptide) const
  {
    if (peptide_.size() != peptide.peptide_.size())
    {
      return false;
    }

    for (Size i = 0; i != size(); ++i)
    {
      if (peptide_[i] != peptide.peptide_[i])
      {
        return false;
      }
      // if AA sequence equal, check if modifications (if available) are equal
      else if (peptide_.at(i)->getModification() != peptide.peptide_.at(i)->getModification())
      {
        return false;
      }
    }

    if (n_term_mod_ != peptide.n_term_mod_)
    {
      return false;
    }

    if (c_term_mod_ != peptide.c_term_mod_)
    {
      return false;
    }

    return true;
  }

  bool AASequence::operator!=(const AASequence& peptide) const
  {
    return !(*this == peptide);
  }

  bool AASequence::empty() const
  {
    return peptide_.empty();
  }

  bool AASequence::isModified() const
  {
    if (n_term_mod_ != nullptr || c_term_mod_ != nullptr)
    {
      return true;
    }

    for (std::vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
    {
      if ((*it)->isModified())
      {
        return true;
      }
    }
    return false;
  }

  std::ostream& operator<<(std::ostream& os, const AASequence& peptide)
  {
    // this is basically the implementation of toString

    // deal with N-terminal modifications first
    if (peptide.n_term_mod_ != nullptr)
    {
      os << peptide.n_term_mod_->toString();
    }

    for (const auto& aa : peptide)
    {
      os << aa.toString();
    }

    // deal with C-terminal modifications
    if (peptide.c_term_mod_ != nullptr)
    {
      os << peptide.c_term_mod_->toString();
    }

    return os;
  }

  String::ConstIterator AASequence::parseModRoundBrackets_(const String::ConstIterator str_it,
                                                           const String& str,
                                                           AASequence& aas,
                                                           const ResidueModification::TermSpecificity& specificity)
  {
    OPENMS_PRECONDITION(*str_it == '(', "Modification must start with '('.");

    String::ConstIterator mod_start = str_it;
    String::ConstIterator mod_end = ++mod_start;
    Size open_brackets = 1;
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    while (mod_end != str.end())
    {
      if (*mod_end == ')')
      {
        --open_brackets;
      }
      else if (*mod_end == '(')
      {
        ++open_brackets;
      }

      if (!open_brackets)
      {
        break;
      }
      ++mod_end;
    }

    // Extract the actual modification as string
    std::string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
          "Cannot convert string to peptide modification: missing ')'");
    }

    // First search for N or C terminal modifications (start of peptide indicates N-terminal modification as well)
    if (aas.peptide_.empty() || specificity == ResidueModification::N_TERM ||
                                specificity == ResidueModification::PROTEIN_N_TERM)
    {
      // Advance iterator one or two positions (we may or may not have a dot
      // after the closing bracket) to point to the first AA of the peptide.
      String::ConstIterator next_aa = mod_end;
      ++next_aa;
      if (*next_aa == '.') ++next_aa;

      if (specificity == ResidueModification::PROTEIN_N_TERM)
      {
        aas.n_term_mod_ = proteinTerminalResidueHelper(mod_db, 'n', str, mod, String(*next_aa));
        return mod_end;
      }
      else
      {
        //TODO why are we allowing Protein Term here?
        aas.n_term_mod_ = terminalResidueHelper(mod_db, 'n', true, str, mod, String(*next_aa));
        return mod_end;
      }
    }

    // get one letter code of unmodified version
    const String& res = aas.peptide_.back()->getOneLetterCode();
    if (specificity == ResidueModification::PROTEIN_C_TERM)
    {
      aas.c_term_mod_ = proteinTerminalResidueHelper(mod_db, 'c', str, mod, res);
      return mod_end;
    }
    else if (specificity == ResidueModification::C_TERM)
    {
      //TODO why are we allowing Protein Term here?
      aas.c_term_mod_ = terminalResidueHelper(mod_db, 'c', true, str, mod, res);
      return mod_end;
    }

    try
    {
      const Residue* internal = ResidueDB::getInstance()->getModifiedResidue(aas.peptide_.back(), mod);
      aas.peptide_.back() = internal;
    }
    catch(...)  // no internal mod for this residue
    {
      // TODO: get rid of this code path, its deprecated and is only a hack for
      // C/N-terminal modifications that don't use the dot notation
      if (std::distance(str_it, str.begin()) == -1)
      {
        // old ambiguous notation: Modification might be at first amino acid or at N-terminus
        aas.n_term_mod_ = terminalResidueHelper(mod_db, 'n', true, str, mod, res);
      }
      else if (std::distance(mod_end, str.end()) == 1) // potentially a C-terminal mod without explicitly declaring it using dot notation?
      {
        // old ambiguous notation: Modification might be at last amino acid or at C-terminus
        aas.c_term_mod_ = terminalResidueHelper(mod_db, 'c', true, str, mod, res);
      }
      else
      {
        // neither internal nor terminal modification matches to our database
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Cannot convert string to peptide modification. No modification matches in our database.", mod);
      }
    }
    return mod_end;
  }

  String::ConstIterator AASequence::parseModSquareBrackets_(const String::ConstIterator str_it,
                                                            const String& str,
                                                            AASequence& aas,
                                                            const ResidueModification::TermSpecificity& specificity)
  {
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");

    String::ConstIterator mod_start = str_it;
    String::ConstIterator mod_end = ++mod_start;
    while ((mod_end != str.end()) && (*mod_end != ']')) {++mod_end;}

    // Extract the actual modification as string
    String mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
          "Cannot convert string to peptide modification: missing ']'");
    }

    double mass = mod.toDouble();
    size_t decimal_pos = mod.find('.');
    bool integer_mass = (decimal_pos == std::string::npos);
    double tolerance = 0.5; // for integer mass values
    if (!integer_mass) // float mass values -> adapt tolerance to decimal precision
    {
      size_t n_decimals = mod.size() - decimal_pos - 2;
      tolerance = std::pow(10.0, -int(n_decimals));
    }
    bool delta_mass = (mod[0] == '+') || (mod[0] == '-');
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    const Residue* residue = nullptr;
    const ResidueModification* residue_mod = nullptr;

    // handle N-term modification
    if (specificity == ResidueModification::N_TERM)
    {
      // Advance iterator one or two positions (we may or may not have a dot
      // after the closing bracket) to point to the first AA of the peptide.
      String::ConstIterator next_aa = mod_end;
      ++next_aa;
      if (*next_aa == '.') ++next_aa;
      std::vector<String> term_mods;
      if (!integer_mass) // for non-integer mass we just pick the closest inside the tolerance
      {
        if (delta_mass) // N-terminal mod specified by delta mass [+123.4]
        {
          residue_mod = mod_db->getBestModificationByDiffMonoMass(mass, tolerance, String(*next_aa), ResidueModification::N_TERM);
        }
        else // N-terminal mod specified by absolute mass [123.4]
        {
          double mod_mass = mass - Residue::getInternalToNTerm().getMonoWeight(); // here we need to subtract the N-Term mass
          residue_mod = mod_db->getBestModificationByDiffMonoMass(mod_mass, tolerance, String(*next_aa), ResidueModification::N_TERM);
        }
        if (residue_mod != nullptr)
        {
          aas.n_term_mod_ = residue_mod;
          return mod_end;
        }
      }
      else // for integer mass we report on the modification in the tolerance and report which we picked.
      {
        if (delta_mass) // N-terminal mod specified by delta mass [+123.4]
        {
          mod_db->searchModificationsByDiffMonoMass(term_mods, mass, tolerance, String(*next_aa), ResidueModification::N_TERM);
        }
        else // N-terminal mod specified by absolute mass [123.4]
        {
          double mod_mass = mass - Residue::getInternalToNTerm().getMonoWeight(); // here we need to subtract the N-Term mass
          mod_db->searchModificationsByDiffMonoMass(term_mods, mod_mass, tolerance, String(*next_aa), ResidueModification::N_TERM);
        }
        if (!term_mods.empty())
        {
          aas.n_term_mod_ = mod_db->getModification(term_mods[0], String(*next_aa), ResidueModification::N_TERM);
          return mod_end;
        }
      }

      OPENMS_LOG_WARN << "Warning: unknown N-terminal modification '" + mod + "' - adding it to the database" << std::endl;
    }
    else if (specificity == ResidueModification::ANYWHERE) // internal (not exclusively terminal) modification
    {
      residue = aas.peptide_.back();
      if (delta_mass && (residue->getMonoWeight() <= 0.0)) // not allowed
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
            "Using a mass difference to specify a modification on a residue of unknown mass is not supported in '" + \
            residue->getOneLetterCode() + "[" + mod +  "]'");
      }

      if (!delta_mass) // compute delta mass based on residue mass
      {
        // we expect that all masses are relative to the full mass of the
        // residue, not its "internal" mass in the peptide (loss of H2O)!
        mass -= residue->getMonoWeight(Residue::Internal);
        delta_mass = true; // in case we need to create a new residue below
      }

      if (integer_mass) // use first modification that matches approximately
      {
        std::vector<String> res_mods;

        mod_db->searchModificationsByDiffMonoMass(res_mods, mass, tolerance, residue->getOneLetterCode(), ResidueModification::ANYWHERE);
        if (!res_mods.empty())
        {
          aas.peptide_.back() = ResidueDB::getInstance()->getModifiedResidue(residue, res_mods[0]);
          return mod_end;
        }
        else if (aas.size() == 1) // N-terminal mod.?
        {
          std::vector<String> term_mods;
          mod_db->searchModificationsByDiffMonoMass(term_mods, mass, tolerance, residue->getOneLetterCode(), ResidueModification::N_TERM);
          if (!term_mods.empty())
          {
            aas.n_term_mod_ = mod_db->getModification(term_mods[0], residue->getOneLetterCode(), ResidueModification::N_TERM);
            return mod_end;
          }
        }
        else if (std::distance(mod_end, str.end()) == 1) // C-terminal mod.?
        {
          mod_db->searchModificationsByDiffMonoMass(res_mods, mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
          if (!res_mods.empty())
          {
            aas.c_term_mod_ = mod_db->getModification(res_mods[0], residue->getOneLetterCode(), ResidueModification::C_TERM);
            return mod_end;
          }
        }
      }
      else // float mass -> use best-matching modification
      {
        const ResidueModification* res_mod = mod_db->getBestModificationByDiffMonoMass(
          mass, tolerance, residue->getOneLetterCode(),
          ResidueModification::ANYWHERE);

        if (res_mod)
        {
          String id = res_mod->getId();
          if (id.empty()) id = res_mod->getFullId();
          aas.peptide_.back() = ResidueDB::getInstance()->getModifiedResidue(residue, id);
          return mod_end;
        }
        else if (aas.size() == 1) // N-terminal mod.?
        {
          res_mod = mod_db->getBestModificationByDiffMonoMass(mass, tolerance, residue->getOneLetterCode(), ResidueModification::N_TERM);
          if (res_mod)
          {
            aas.n_term_mod_ = res_mod;
            return mod_end;
          }
        }
        else if (std::distance(mod_end, str.end()) == 1) // C-terminal mod.?
        {
          res_mod = mod_db->getBestModificationByDiffMonoMass(mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
          if (res_mod)
          {
            aas.c_term_mod_ = res_mod;
            return mod_end;
          }
        }
      }
      if (residue->getOneLetterCode() != "X") // don't warn for mass tags
      {
        OPENMS_LOG_WARN << "Warning: unknown modification '" + mod + "' of residue '" +
            residue->getOneLetterCode() + "' - adding it to the database" << std::endl;
      }
    }
    else if (specificity == ResidueModification::C_TERM)
    {
      residue = aas.peptide_.back();
      std::vector<String> term_mods;

      if (!integer_mass) // for non-integer mass we just pick the closest inside the tolerance
      {
        if (delta_mass) // C-terminal mod specified by delta mass [+123.4]
        {
          residue_mod = mod_db->getBestModificationByDiffMonoMass(mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
        }
        else // C-terminal mod specified by absolute mass [123.4]
        {
          double mod_mass = mass - Residue::getInternalToCTerm().getMonoWeight(); // here we need to subtract the N-Term mass
          residue_mod = mod_db->getBestModificationByDiffMonoMass(mod_mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
        }
        if (residue_mod != nullptr)
        {
          aas.c_term_mod_ = residue_mod;
          return mod_end;
        }
      }
      else // for integer mass we report on the modification in the tolerance and report which we picked.
      {
        if (delta_mass) // C-terminal mod specified by delta mass [+123]
        {
          mod_db->searchModificationsByDiffMonoMass(term_mods, mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
        }
        else // C-terminal mod specified by absolute mass [123]
        {
          double mod_mass = mass - Residue::getInternalToCTerm().getMonoWeight(); // here we need to subtract the N-Term mass
          mod_db->searchModificationsByDiffMonoMass(term_mods, mod_mass, tolerance, residue->getOneLetterCode(), ResidueModification::C_TERM);
        }
        if (!term_mods.empty())
        {
          aas.c_term_mod_ = mod_db->getModification(term_mods[0], residue->getOneLetterCode(), ResidueModification::C_TERM);
          return mod_end;
        }
      }

      OPENMS_LOG_WARN << "Warning: unknown C-terminal modification '" + mod + "' - adding it to the database" << std::endl;
    }

    // -----------------------------------
    // Dealing with an unknown modification
    // -----------------------------------

    const ResidueModification* new_mod = ResidueModification::createUnknownFromMassString(mod,
                                                                                          mass,
                                                                                          delta_mass,
                                                                                          specificity,
                                                                                          residue);

    // Notes on mass calculation: AASequence::getMonoWeight uses DiffMonoMass
    // for its calculation of C/N-terminal modification mass and it uses
    // getMonoWeight(Residue::Internal) for each Residue. The Residue weight is
    // set when adding a modification using setModification_
    if (specificity == ResidueModification::N_TERM)
    {
      aas.n_term_mod_ = new_mod;
      return mod_end;
    }
    else if (specificity == ResidueModification::C_TERM)
    {
      aas.c_term_mod_ = new_mod;
      return mod_end;
    }
    else
    {
      // Note: this calls setModification_ on a new Residue which changes its
      // weight to the weight of the modification (set above)
      aas.peptide_.back() = ResidueDB::getInstance()->
        getModifiedResidue(residue, new_mod->getFullId());
      return mod_end;
    }
  }

  void AASequence::parseString_(const String& pep, AASequence& aas,
                                bool permissive)
  {
    // Reserving space, populate it and then shrink again (since we probably
    // over-allocate due to modifications). This substantially speeds up the
    // function for unmodified sequences (3x speedup).
    aas.peptide_.clear();

    String peptide(pep);
    peptide.trim();
    aas.peptide_.reserve(peptide.size());

    if (peptide.empty()) return;

    // remove optional n and c at start and end of string
    if (peptide[0] == 'n')
    {
      peptide.erase(0,1);
    }
    if (peptide.empty())
    {
      return;
    }
    if (peptide.back() == 'c')
    {
      peptide.pop_back();
    }

    if (peptide.empty()) return;

    // detect if this is the new dot notation containing dots for termini and
    // track if last char denoted a terminus
    bool dot_terminal(false), dot_notation(false);

    static ResidueDB* rdb = ResidueDB::getInstance();

    for (String::ConstIterator str_it = peptide.begin();
         str_it != peptide.end(); ++str_it)
    {
      // skip (optional) terminal delimiters, but remember that last character was a terminal one
      if (*str_it == '.')
      {
        dot_notation = true;
        dot_terminal = true;
        continue;
      }

      // 1. default case: add unmodified, standard residue
      const Residue* r = rdb->getResidue(*str_it); // "isalpha" check not needed
      if (r)
      {
        dot_terminal = false; // since we found an AA, we are not at a terminal position any more
        aas.peptide_.push_back(r);
        continue;
      }

      // 2. modification:
      //   determine specificity:
      //     - at termini we first assume we are dealing with a N- or C-terminal modifications
      //       and fall back to (internal) modifications if there is none in our DB
      //     - otherwise we can be sure we are dealing with an internal modification
      ResidueModification::TermSpecificity specificity = ResidueModification::ANYWHERE;

      if (str_it == peptide.begin() || (dot_notation && dot_terminal && aas.peptide_.empty()) )
      {
        specificity = ResidueModification::N_TERM;
      }
      else if (*str_it == 'c')
      {
        // note that still c[...] type substring remains as only single c have been erased before
        // skip 'c', record that we are dealing with a C-terminal
        ++str_it;
        specificity = ResidueModification::C_TERM;
      }
      else if (dot_notation && dot_terminal && !aas.peptide_.empty())
      {
        specificity = ResidueModification::C_TERM;
      }

      if (*str_it == '(')
      {
        str_it = parseModRoundBrackets_(str_it, peptide, aas, specificity);
      }
      else if (*str_it == '[')
      {
        str_it = parseModSquareBrackets_(str_it, peptide, aas, specificity);
      }
      else
      {
        if (permissive && ((*str_it == '*') || (*str_it == '#') ||
                           (*str_it == '+')))
        { // stop codons
          aas.peptide_.push_back(rdb->getResidue('X'));
        }
        else if (permissive && (*str_it == ' '))
        { // skip, i.e. do nothing here
        }
        else
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, peptide,
              "Cannot convert string to amino acid sequence: unexpected character '" + String(*str_it) + "'");
        }
      }

      dot_terminal = false; // previous char was no dot
    }

    // We do NOT deal with a single, unmodified X residue here,
    // since the user might just want to represent the sequence (including modifications on other AA's),
    // e.g. when digesting a peptide
    // We check for 'weightless' X in places where a mass is needed, e.g. during getMonoMass()

    aas.peptide_.shrink_to_fit();
  }

  void AASequence::getAAFrequencies(std::map<String, Size>& frequency_table) const
  {
    frequency_table.clear();
    for (std::vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
    {
      ++frequency_table[(*it)->getOneLetterCode()];
    }
  }

  void AASequence::setModification(Size index, const String& modification)
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, peptide_.size());
    }

    if (!modification.empty())
    {
      peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], modification);
    }
    else // remove modification
    {
      peptide_[index] = ResidueDB::getInstance()->getResidue(peptide_[index]->getOneLetterCode());
    }
  }

  void AASequence::setModification(Size index, const ResidueModification* modindb)
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, peptide_.size());
    }
    peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], modindb);
  }

  void AASequence::setModification(Size index, const ResidueModification& modification)
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, peptide_.size());
    }
    const auto& db = ModificationsDB::getInstance();
    const ResidueModification* modindb = db->searchModification(modification);
    if (modindb == nullptr)
    {
      modindb = db->addNewModification_(modification);
    }
    peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], modindb);
  }

  void AASequence::setModificationByDiffMonoMass(Size index, double diffMonoMass)
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, peptide_.size());
    }

    ModificationsDB* mod_db = ModificationsDB::getInstance();
    bool multimatch = false;
    // quickly check for user-defined modification added by createUnknownFromMassString (e.g. M[+12321])
    String diffMonoMassStr = ResidueModification::getDiffMonoMassWithBracket(diffMonoMass);
    const ResidueModification* mod = mod_db->searchModificationsFast(peptide_[index]->getOneLetterCode() + diffMonoMassStr, multimatch);
    const double tol = 0.002;
    if (mod == nullptr)
    {
      mod = mod_db->getBestModificationByDiffMonoMass(diffMonoMass, tol, peptide_[index]->getOneLetterCode(), ResidueModification::ANYWHERE);
    }
    if (mod == nullptr)
    {
      OPENMS_LOG_WARN << "Modification with monoisotopic mass diff. of " << diffMonoMassStr << " not found in databases with tolerance " << tol << ". Adding unknown modification." << std::endl;
      mod = ResidueModification::createUnknownFromMassString(String(diffMonoMass),
                                                                        diffMonoMass,
                                                                        true,
                                                                        ResidueModification::ANYWHERE,
                                                                        peptide_[index]);
    }
    peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], mod);
  }

  void AASequence::setModification(Size index, const Residue* modification)
  {
    peptide_[index] = modification;
  }

  void AASequence::setNTerminalModification(const String& modification)
  {
    if (modification.empty())
    {
      n_term_mod_ = nullptr;
      return;
    }

    String residue = "";
    if (modification.size() > 3 && modification.hasSuffix(")"))
    {
      char last_char_no_parentheses = modification[modification.length()-2];
      if (isupper(last_char_no_parentheses))
      {
        residue = last_char_no_parentheses;
      }
    }

    // For strings in our most common UniMod format
    if (!modification.hasSubstring("Protein N-term"))
    {
      // since this method is called setNTerminalModification without further specification
      // we have to look for both general terminus and Protein terminus.
      // For backwards compatibility we look for the general terminus first
      // We have to use try-catch since getModification unfortunately throws Exceptions.
      try
      {
        n_term_mod_ = ModificationsDB::getInstance()
            ->getModification(modification, residue, ResidueModification::N_TERM);
        return; // we found a mod. return
      }
      catch (...) {}
    }
    n_term_mod_ = ModificationsDB::getInstance()->getModification(modification, residue, ResidueModification::PROTEIN_N_TERM);
  }

  void AASequence::setCTerminalModification(const String& modification)
  {
    if (modification.empty())
    {
      c_term_mod_ = nullptr;
      return;
    }

    String residue = "";
    if (modification.size() > 3 && modification.hasSuffix(")"))
    {
      char last_char_no_parentheses = modification[modification.length()-2];
      if (isupper(last_char_no_parentheses))
      {
        residue = last_char_no_parentheses;
      }
    }

    // For strings in our most common UniMod format
    if (!modification.hasSubstring("Protein C-term"))
    {
      // since this method is called setCTerminalModification without further specification
      // we have to look for both general terminus and Protein terminus.
      // For backwards compatibility we look for the general terminus first
      // We have to use try-catch since getModification unfortunately throws Exceptions.
      try
      {
        c_term_mod_ = ModificationsDB::getInstance()->getModification(modification, residue, ResidueModification::C_TERM);
        return; // we found a mod. return
      }
      catch (...) {}
    }
    c_term_mod_ = ModificationsDB::getInstance()->getModification(modification, residue, ResidueModification::PROTEIN_C_TERM);
  }

  void AASequence::setCTerminalModification(const ResidueModification& mod)
  {
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    //TODO think again. Most functions here or in ModificationsDB only check for fullID
    c_term_mod_ = mod_db->searchModification(mod);
    if (c_term_mod_ == nullptr)
    {
      c_term_mod_ = mod_db->addNewModification_(mod);
    }
  }

  void AASequence::setNTerminalModification(const ResidueModification& mod)
  {
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    //TODO think again. Most functions here or in ModificationsDB only check for fullID
    n_term_mod_ = mod_db->searchModification(mod);
    if (n_term_mod_ == nullptr)
    {
      n_term_mod_ = mod_db->addNewModification_(mod);
    }
  }

  void AASequence::setCTerminalModification(const ResidueModification* modification)
  {
    c_term_mod_ = modification;
  }

  void AASequence::setNTerminalModification(const ResidueModification* modification)
  {
    n_term_mod_ = modification;
  }

  void AASequence::setCTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term)
  {
    // since this method is called setCTerminalModification without further specification
    // we have to look for both general terminus and Protein terminus.
    // For backwards compatibility we look for the general terminus first
    ResidueModification::TermSpecificity term = protein_term ? ResidueModification::PROTEIN_C_TERM : ResidueModification::C_TERM;
    double tol = 0.002;
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    bool multimatch = false;
    // quickly check for user-defined modification added by createUnknownFromMassString (e.g. M[+12321])
    String diffMonoMassStr = ResidueModification::getDiffMonoMassWithBracket(diffMonoMass);
    // TODO make a distinction in the FullID about protein vs peptide term??
    const ResidueModification* n_term_mod_ = mod_db->searchModificationsFast(".c"+diffMonoMassStr, multimatch);
    String residue = "";
    if (n_term_mod_ == nullptr)
    {
        n_term_mod_ = ModificationsDB::getInstance()
          ->getBestModificationByDiffMonoMass(diffMonoMass, tol, residue, term);
    }
    if (n_term_mod_ == nullptr)
    {

      OPENMS_LOG_WARN << "Modification with monoisotopic mass diff. of " << diffMonoMassStr << " not found in databases with tolerance " << tol << ". Adding unknown modification." << std::endl;
      n_term_mod_ = ResidueModification::createUnknownFromMassString(String(diffMonoMass),
                                                                        diffMonoMass,
                                                                        true,
                                                                        term);
    }
  }

  void AASequence::setNTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term)
  {
    // since this method is called setNTerminalModification without further specification
    // we have to look for both general terminus and Protein terminus.
    // For backwards compatibility we look for the general terminus first
    ResidueModification::TermSpecificity term = protein_term ? ResidueModification::PROTEIN_N_TERM : ResidueModification::N_TERM;
    double tol = 0.002;
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    bool multimatch = false;
    // quickly check for user-defined modification added by createUnknownFromMassString (e.g. M[+12321])
    String diffMonoMassStr = ResidueModification::getDiffMonoMassWithBracket(diffMonoMass);
    // TODO make a distinction in the FullID about protein vs peptide term??
    const ResidueModification* n_term_mod_ = mod_db->searchModificationsFast(".n"+diffMonoMassStr, multimatch);
    String residue = "";
    if (n_term_mod_ == nullptr)
    {
        n_term_mod_ = ModificationsDB::getInstance()
          ->getBestModificationByDiffMonoMass(diffMonoMass, tol, residue, term);
    }
    if (n_term_mod_ == nullptr)
    {

      OPENMS_LOG_WARN << "Modification with monoisotopic mass diff. of " << diffMonoMassStr << " not found in databases with tolerance " << tol << ". Adding unknown modification." << std::endl;
      n_term_mod_ = ResidueModification::createUnknownFromMassString(String(diffMonoMass),
                                                                        diffMonoMass,
                                                                        true,
                                                                        term);
    }
  }

  const String& AASequence::getNTerminalModificationName() const
  {
    if (n_term_mod_ == nullptr) return String::EMPTY;
    return n_term_mod_->getId();
  }

  const ResidueModification* AASequence::getNTerminalModification() const
  {
    return n_term_mod_;
  }

  const ResidueModification* AASequence::getCTerminalModification() const
  {
    return c_term_mod_;
  }

  const String& AASequence::getCTerminalModificationName() const
  {
    if (c_term_mod_ == nullptr) return String::EMPTY;
    return c_term_mod_->getId();
  }

  bool AASequence::hasNTerminalModification() const
  {
    return n_term_mod_ != nullptr;
  }

  bool AASequence::hasCTerminalModification() const
  {
    return c_term_mod_ != nullptr;
  }

  AASequence AASequence::fromString(const String& s, bool permissive)
  {
    AASequence aas;
    parseString_(s, aas, permissive);
    return aas;
  }

  AASequence AASequence::fromString(const char* s, bool permissive)
  {
    AASequence aas;
    parseString_(String(s), aas, permissive);
    return aas;
  }

}
