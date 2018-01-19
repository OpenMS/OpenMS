// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <string>
#include <algorithm>
#include <cmath> // for "pow"
#include <iterator> // for "distance"
#include <sstream>

using namespace std;

namespace OpenMS
{
  AASequence::AASequence() :
    n_term_mod_(nullptr),
    c_term_mod_(nullptr)
  {
  }

  AASequence::AASequence(const AASequence& rhs) :
    peptide_(rhs.peptide_),
    n_term_mod_(rhs.n_term_mod_),
    c_term_mod_(rhs.c_term_mod_)
  {
  }

  AASequence::~AASequence()
  {
  }

  AASequence& AASequence::operator=(const AASequence& rhs)
  {
    if (this != &rhs)
    {
      peptide_ = rhs.peptide_;
      n_term_mod_ = rhs.n_term_mod_;
      c_term_mod_ = rhs.c_term_mod_;
    }
    return *this;
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
      if (mod.getUniModRecordId() > -1)
      {
        bs += ".(" + mod.getUniModAccession() + ")";
      }
      else
      {
        bs += ".[" + String(mod.getDiffMonoMass()) + "]";
      }
    }

    for (Size i = 0; i != seq.size(); ++i)
    {
      const Residue& r = seq[i];
      const String aa = r.getOneLetterCode();
      if (r.isModified())
      {
        const ResidueModification& mod = *(r.getModification());
        if (mod.getUniModRecordId() > -1)
        {
          bs += aa + "(" + mod.getUniModAccession() + ")";
        }
        else
        {
          bs += aa + "[" + String(r.getMonoWeight(Residue::Internal)) + "]";
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
      if (mod.getUniModRecordId() > -1)
      {
        bs += ".(" + mod.getUniModAccession() + ")";
      }
      else
      {
        bs += ".[" + String(mod.getDiffMonoMass()) + "]";
      }
    }
    return bs;
  }

  String AASequence::toBracketString(bool integer_mass, const vector<String> & fixed_modifications) const
  {
    const AASequence & seq = *this;

    String bs;
    if (seq.empty()) return bs;

    if (seq.hasNTerminalModification())
    {
      const ResidueModification& mod = *(seq.getNTerminalModification());
      const String & nterm_mod_name = mod.getFullId(); // e.g. "Acetyl (N-term)"

      // only add to string if not a fixed modification
      if (std::find(fixed_modifications.begin(), fixed_modifications.end(), nterm_mod_name) == fixed_modifications.end())
      {
        double nominal_mass = Residue::getInternalToNTerm().getMonoWeight() + mod.getDiffMonoMass();
        if (mod.isUserDefined()) nominal_mass = mod.getDiffMonoMass(); // just get input mass
        if (integer_mass)
        {
          bs += "n[" + String(static_cast<int>(nominal_mass)) + "]";
        }
        else
        {
          bs += "n[" + String(nominal_mass) + "]";
        }
      }
    }

    for (Size i = 0; i != seq.size(); ++i)
    {
      const Residue& r = seq[i];
      const String aa = r.getOneLetterCode() != "" ? r.getOneLetterCode() : "X";
      if (r.isModified())
      {
        const ResidueModification& mod = *(r.getModification());

        const String & mod_name = mod.getFullId();
        if (std::find(fixed_modifications.begin(), fixed_modifications.end(), mod_name) == fixed_modifications.end())
        {
          if (integer_mass)
          {
            const double residue_mono_mass = r.getMonoWeight(Residue::Internal);
            bs += aa + "[" + static_cast<int>(residue_mono_mass) + "]"; 
          }
          else
          {
            const double residue_mono_mass = r.getMonoWeight(Residue::Internal);
            bs += aa + "[" + residue_mono_mass + "]"; 
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
        double nominal_mass = Residue::getInternalToCTerm().getMonoWeight() + mod.getDiffMonoMass();
        if (mod.isUserDefined()) nominal_mass = mod.getDiffMonoMass(); // just get input mass
        if (integer_mass)
        {
          bs += "c[" + String(static_cast<int>(nominal_mass)) + "]";
        }
        else
        {
          bs += "c[" + String(nominal_mass) + "]";
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
    OPENMS_PRECONDITION( !this->has( *ResidueDB::getInstance()->getResidue("X") ), "Cannot get formula of sequence with unknown AA with unknown mass.");

    if (peptide_.size() >= 1)
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

      for (Size i = 0; i != peptide_.size(); ++i)
      {
        ef += peptide_[i]->getFormula(Residue::Internal);
      }

          // add the missing formula part
      switch (type)
      {
        case Residue::Full:
          return ef + Residue::getInternalToFull();

        case Residue::Internal:
          return ef;

        case Residue::NTerminal:
          return ef + Residue::getInternalToNTerm();

        case Residue::CTerminal:
          return ef + Residue::getInternalToCTerm();

        case Residue::AIon:
          return ef + Residue::getInternalToAIon();

        case Residue::BIon:
          return ef + Residue::getInternalToBIon();

        case Residue::CIon:
          return ef + Residue::getInternalToCIon();

        case Residue::XIon:
          return ef + Residue::getInternalToXIon();

        case Residue::YIon:
          return ef + Residue::getInternalToYIon();

        case Residue::ZIon:
          return ef + Residue::getInternalToZIon();

        default:
          LOG_ERROR << "AASequence::getFormula: unknown ResidueType" << std::endl;
      }

      return ef;
    }
    else
    {
      LOG_ERROR << "AASequence::getFormula: Formula for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
      return EmpiricalFormula("");
    }
  }

  double AASequence::getAverageWeight(Residue::ResidueType type, Int charge) const
  {
    OPENMS_PRECONDITION( !this->has( *ResidueDB::getInstance()->getResidue("X") ), "Cannot get weight of sequence with unknown AA with unknown mass.");

    // check whether tags are present
    double tag_offset(0);
    for (ConstIterator it = this->begin(); it != this->end(); ++it)
    {
      if (it->getOneLetterCode() == "")
      {
        tag_offset += it->getAverageWeight(Residue::Internal);
      }
    }
    // TODO inefficient, if averageWeight is already set in the Residue
    return tag_offset + getFormula(type, charge).getAverageWeight();
  }

  double AASequence::getMonoWeight(Residue::ResidueType type, Int charge) const
  {
    OPENMS_PRECONDITION( !this->has( *ResidueDB::getInstance()->getResidue("X") ), "Cannot get weight of sequence with unknown AA with unknown mass.");

    if (peptide_.size() >= 1)
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

      for (ConstIterator it = this->begin(); it != this->end(); ++it)
      {
        // standard internal residue including named modifications
        mono_weight += it->getMonoWeight(Residue::Internal);
      }

      // add the missing formula part
      switch (type)
      {
        case Residue::Full:
          return mono_weight + Residue::getInternalToFull().getMonoWeight();

        case Residue::Internal:
          return mono_weight;

        case Residue::NTerminal:
          return mono_weight + Residue::getInternalToNTerm().getMonoWeight();

        case Residue::CTerminal:
          return mono_weight + Residue::getInternalToCTerm().getMonoWeight();

        case Residue::AIon:
          return mono_weight + Residue::getInternalToAIon().getMonoWeight();

        case Residue::BIon:
          return mono_weight + Residue::getInternalToBIon().getMonoWeight();

        case Residue::CIon:
          return mono_weight + Residue::getInternalToCIon().getMonoWeight();

        case Residue::XIon:
          return mono_weight + Residue::getInternalToXIon().getMonoWeight();

        case Residue::YIon:
          return mono_weight + Residue::getInternalToYIon().getMonoWeight();

        case Residue::ZIon:
          return mono_weight + Residue::getInternalToZIon().getMonoWeight();

        default:
          LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
      }

      return mono_weight;
    }
    else
    {
      LOG_ERROR << "AASequence::getMonoWeight: Mass for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
      return 0.0;
    }
}




/*void AASequence::getNeutralLosses(Map<const EmpiricalFormula, UInt) const
  {
      // the following losses are from the Zhang paper (AC, 76, 14, 2004)
      // charge directed*/
/*
  static const EmpiricalFormula R_44("NH2CHNH");
  static const EmpiricalFormula R_59("CN3H5"); // guanidium
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
      seq.n_term_mod_ = n_term_mod_;
    if (index + num == this->size())
      seq.c_term_mod_ = c_term_mod_;

    seq.peptide_.insert(seq.peptide_.end(), peptide_.begin() + index, peptide_.begin() + index + num);

    return seq;
  }

  bool AASequence::has(const Residue& residue) const
  {
    for (Size i = 0; i != peptide_.size(); ++i)
    {
      if (*peptide_[i] == residue)
      {
        return true;
      }
    }
    return false;
  }

  bool AASequence::hasSubsequence(const AASequence& sequence) const
  {
    if (sequence.empty())
    {
      return true;
    }
    else
    {
      if (sequence.size() <= peptide_.size())
      {
        for (Size i = 0; i != peptide_.size(); ++i)
        {
          if (peptide_[i] == sequence.peptide_[0])
          {
            Size j = 0;
            for (; j + i != peptide_.size() && j != sequence.peptide_.size(); ++j)
            {
              if (peptide_[j + i] == sequence.peptide_[j])
              {
                if (j == sequence.peptide_.size() - 1)
                {
                  return true;
                }
              }
              else
              {
                break;
              }
            }
          }
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
      return false;

    if (sequence.size() == peptide_.size() && sequence.c_term_mod_ != c_term_mod_)
      return false;

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
      return false;

    if (sequence.size() == peptide_.size() && sequence.n_term_mod_ != n_term_mod_)
      return false;

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
    return size() == 0;
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
      if (peptide.n_term_mod_->isUserDefined())
      {
        os << peptide.n_term_mod_->getFullId();
      }
      else
      {
        os << ".(" << peptide.n_term_mod_->getId() << ")";
      }
    }

    for (Size i = 0; i != peptide.size(); ++i)
    {
      if (peptide.peptide_[i]->isModified())
      {
        if (peptide.peptide_[i]->getOneLetterCode() != "")
        {
          os << peptide.peptide_[i]->getOneLetterCode();
        }
        else
        {
          OPENMS_PRECONDITION(false, "Something went wrong, this should never happen"); // we should not end up here any more
          os << "[" << precisionWrapper(peptide.peptide_[i]->getMonoWeight()) << "]";
        }
        const String& id = peptide.peptide_[i]->getModificationName();
        if (peptide.peptide_[i]->getModification()->isUserDefined())
        {
          // user-defined modification
          os << peptide.peptide_[i]->getModification()->getFullId();
        }
        else if (id != "")
        {
          os << "(" << id << ")";
        }
        else
        {
          // what is this syntax? we have no idea, we probably should never be here
          OPENMS_PRECONDITION(false, "Something went wrong, this should never happen"); // we should not end up here any more
          os << "([" << precisionWrapper(peptide.peptide_[i]->getModification()->getDiffMonoMass()) << "])";
        }
      }
      else
      {
        if (peptide.peptide_[i]->getOneLetterCode() != "")
        {
          os << peptide.peptide_[i]->getOneLetterCode();
        }
        else
        {
          if (peptide.peptide_[i]->getShortName() != "")
          {
            os << peptide.peptide_[i]->getShortName();
          }
          else
          {
            OPENMS_PRECONDITION(false, "Something went wrong, this should never happen"); // we should not end up here any more
            os << "[" << precisionWrapper(peptide.peptide_[i]->getMonoWeight()) << "]";
          }
        }
      }
    }
    
    // deal with C-terminal modifications
    if (peptide.c_term_mod_ != nullptr)
    {
      if (peptide.c_term_mod_->isUserDefined())
      {
        os << peptide.c_term_mod_->getFullId();
      }
      else
      {
        os << ".(" << peptide.c_term_mod_->getId() << ")";
      }
    }

    return os;
  }

  String::ConstIterator AASequence::parseModRoundBrackets_(
    const String::ConstIterator str_it, const String& str, AASequence& aas, bool dot_notation, bool dot_terminal)
  {
    OPENMS_PRECONDITION(*str_it == '(', "Modification must start with '('.");
    String::ConstIterator mod_start = str_it;
    String::ConstIterator mod_end = ++mod_start;
    Size open_brackets = 1;
    while (mod_end != str.end())
    {
      if (*mod_end == ')') --open_brackets;
      else if (*mod_end == '(') ++open_brackets;
      if (!open_brackets) break;
      ++mod_end;
    }
    std::string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
          "Cannot convert string to peptide modification: missing ')'");
    }
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    if (aas.peptide_.empty()) // start of peptide -> N-terminal modification?
    {
      aas.n_term_mod_ = &(mod_db->getModification(mod, "",
                                                    ResidueModification::N_TERM));
      return mod_end;
    }

    const String& res = aas.peptide_.back()->getOneLetterCode();
    try
    {
      aas.peptide_.back() = ResidueDB::getInstance()->
        getModifiedResidue(aas.peptide_.back(), mod);
    }
    catch (Exception::InvalidValue) // no such mod for this residue
    {
      if (std::distance(mod_end, str.end()) == 1) // C-terminal mod.?
      {
        if (dot_notation)
        {
          if (dot_terminal)
          {
            const ResidueModification* term_mod =
              &(mod_db->getModification(mod, res, ResidueModification::C_TERM));
            aas.c_term_mod_ = term_mod;
          }
        }
        else // old ambiguous notation: Modification might be at last amino acid or at C-terminus
        {
          try
          {
            // this might throw ElementNotFound, but so be it:
            const ResidueModification* term_mod =
              &(mod_db->getModification(mod, res, ResidueModification::C_TERM));
            aas.c_term_mod_ = term_mod;
          }
          catch (Exception::ElementNotFound& /* e */)
          { // just do nothing, the mod is presumably a non-terminal one
          }
        }
      }
      else throw; // re-throw the InvalidValue
    }

    return mod_end;
  }

  String::ConstIterator AASequence::parseModSquareBrackets_(
    const String::ConstIterator str_it, const String& str, AASequence& aas,
    const ResidueModification::TermSpecificity& specificity)
  {
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
    String::ConstIterator mod_start = str_it;
    String::ConstIterator mod_end = ++mod_start;
    while ((mod_end != str.end()) && (*mod_end != ']')) ++mod_end;
    std::string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
          "Cannot convert string to peptide modification: missing ']'");
    }

    double mass = String(mod).toDouble();
    size_t decimal_pos = mod.find('.');
    bool integer_mass = decimal_pos == std::string::npos;
    double tolerance = 0.5; // for integer mass values
    if (!integer_mass) // float mass values -> adapt tolerance to decimal precision
    {
      size_t n_decimals = mod.size() - decimal_pos - 2;
      tolerance = std::pow(10.0, -int(n_decimals));
    }
    bool delta_mass = (mod[0] == '+') || (mod[0] == '-');
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    const Residue* residue = nullptr;

    // handle N-term modification
    if (specificity == ResidueModification::N_TERM) 
    {
      if (delta_mass) // N-terminal mod specified by delta mass [+123.4]
      {
        std::vector<String> term_mods;
        mod_db->searchModificationsByDiffMonoMass(term_mods, mass, tolerance, "",
                                                  ResidueModification::N_TERM);
        if (!term_mods.empty())
        {
          aas.n_term_mod_ = &(mod_db->getModification(
                                term_mods[0], "", ResidueModification::N_TERM));
          return mod_end;
        }
        LOG_WARN << "Warning: unknown N-terminal modification '" + mod + "' - adding it to the database" << std::endl;
      }
      else // N-terminal mod specified by absolute mass [123.4]
      {
        double mod_mass = mass - Residue::getInternalToNTerm().getMonoWeight(); // here we need to subtract the N-Term mass
        std::vector<String> term_mods;
        mod_db->searchModificationsByDiffMonoMass(term_mods, mod_mass, tolerance, "",
                                                ResidueModification::N_TERM);
        if (!term_mods.empty())
        {
          aas.n_term_mod_ = &(mod_db->getModification(
                                term_mods[0], "", ResidueModification::N_TERM));
          return mod_end;
        }
        LOG_WARN << "Warning: unknown N-terminal modification '" + mod + "' - adding it to the database" << std::endl;
      }
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

        mod_db->searchModificationsByDiffMonoMass(
          res_mods, mass, tolerance, residue->getOneLetterCode(),
          ResidueModification::ANYWHERE);
        if (!res_mods.empty())
        {
          aas.peptide_.back() =
            ResidueDB::getInstance()->getModifiedResidue(residue, res_mods[0]);
          return mod_end;
        }
        else if (std::distance(mod_end, str.end()) == 1) // C-terminal mod.?
        {
          mod_db->searchModificationsByDiffMonoMass(
            res_mods, mass, tolerance, residue->getOneLetterCode(),
            ResidueModification::C_TERM);
          if (!res_mods.empty())
          {
            aas.c_term_mod_ =
              &(mod_db->getModification(res_mods[0],
                                        residue->getOneLetterCode(),
                                        ResidueModification::C_TERM));
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
          aas.peptide_.back() = ResidueDB::getInstance()->
            getModifiedResidue(residue, id);
          return mod_end;
        }
        else if (std::distance(mod_end, str.end()) == 1) // C-terminal mod.?
        {
          res_mod = mod_db->getBestModificationByDiffMonoMass(
            mass, tolerance, residue->getOneLetterCode(),
            ResidueModification::C_TERM);

          if (res_mod)
          {
            aas.c_term_mod_ = res_mod;
            return mod_end;
          }
        }
      }
      LOG_WARN << "Warning: unknown modification '" + mod + "' of residue '" +
        residue->getOneLetterCode() + "' - adding it to the database" << std::endl;
    }
    else if (specificity == ResidueModification::C_TERM)
    {
      if (delta_mass) // C-terminal mod specified by delta mass [+123.4]
      {
        std::vector<String> term_mods;
        mod_db->searchModificationsByDiffMonoMass(term_mods, mass, tolerance, "",
                                                  ResidueModification::C_TERM);
        if (!term_mods.empty())
        {
          aas.c_term_mod_ = &(mod_db->getModification(
                                term_mods[0], "", ResidueModification::C_TERM));
          return mod_end;
        }
        LOG_WARN << "Warning: unknown C-terminal modification '" + mod + "' - adding it to the database" << std::endl;
      }
      else // C-terminal mod specified by absolute mass [123.4]
      {
        double mod_mass = mass - Residue::getInternalToCTerm().getMonoWeight(); // here we need to subtract the C-Term mass
        std::vector<String> term_mods;
        mod_db->searchModificationsByDiffMonoMass(term_mods, mod_mass, tolerance, "",
                                                ResidueModification::C_TERM);
        if (!term_mods.empty())
        {
          aas.c_term_mod_ = &(mod_db->getModification(
                                term_mods[0], "", ResidueModification::C_TERM));
          return mod_end;
        }
        LOG_WARN << "Warning: unknown C-terminal modification '" + mod + "' - adding it to the database" << std::endl;
      }
    }

    // ----------------------------------- 
    // Dealing with an unknown modification
    // ----------------------------------- 

    // Notes on mass calculation: AASequence::getMonoWeight uses DiffMonoMass
    // for its calculation of C/N-terminal modification mass and it uses
    // getMonoWeight(Residue::Internal) for each Residue. The Residue weight is
    // set when adding a modification using setModification_
    if (specificity == ResidueModification::N_TERM) 
    {
      String residue_name = ".[" + mod + "]";

      // Check if it already exists, if not create new modification, transfer
      // ownership to ModDB
      if (!mod_db->has(residue_name)) 
      {
        ResidueModification * new_mod = new ResidueModification();
        new_mod->setFullId(residue_name); // setting FullId but not Id makes it a user-defined mod
        new_mod->setDiffMonoMass(mass);
        new_mod->setTermSpecificity(ResidueModification::N_TERM);
        // new_mod->setMonoMass(mass);
        // new_mod->setAverageMass(mass);
        mod_db->addModification(new_mod);
        aas.n_term_mod_ = new_mod;
      }
      else
      {
        Size mod_idx = mod_db->findModificationIndex(residue_name);
        aas.n_term_mod_ = &mod_db->getModification(mod_idx);
      }
      return mod_end;
    }
    else if (specificity == ResidueModification::C_TERM)
    {
      String residue_name = ".[" + mod + "]";

      // Check if it already exists, if not create new modification, transfer
      // ownership to ModDB
      if (!mod_db->has(residue_name)) 
      {
        ResidueModification * new_mod = new ResidueModification();
        new_mod->setFullId(residue_name); // setting FullId but not Id makes it a user-defined mod
        new_mod->setDiffMonoMass(mass);
        new_mod->setTermSpecificity(ResidueModification::C_TERM);
        // new_mod->setMonoMass(mass);
        // new_mod->setAverageMass(mass);
        mod_db->addModification(new_mod);
        aas.c_term_mod_ = new_mod;
      }
      else
      {
        Size mod_idx = mod_db->findModificationIndex(residue_name);
        aas.c_term_mod_ = &mod_db->getModification(mod_idx);
      }
      return mod_end;
    }
    else
    {
      String residue_name = aas.peptide_.back()->getOneLetterCode() + "[" + mod + "]"; // e.g. N[12345.6]
      String modification_name = "[" + mod + "]";

      if (!mod_db->has(modification_name)) 
      {
        // create new modification
        ResidueModification * new_mod = new ResidueModification();
        new_mod->setFullId(modification_name); // setting FullId but not Id makes it a user-defined mod

        // We cannot set origin if we want to use the same modification name
        // also at other AA (and since we have no information here, it is safer
        // to assume that this may happen).
        // new_mod->setOrigin(aas.peptide_.back()->getOneLetterCode()[0]);

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

        mod_db->addModification(new_mod);
      }

      // now use the new modification
      Size mod_idx = mod_db->findModificationIndex(modification_name);
      const ResidueModification* res_mod = &mod_db->getModification(mod_idx);

      // Note: this calls setModification_ on a new Residue which changes its
      // weight to the weight of the modification (set above)
      aas.peptide_.back() = ResidueDB::getInstance()->
        getModifiedResidue(residue, res_mod->getFullId());
      return mod_end;
    }
  }

  void AASequence::parseString_(const String& pep, AASequence& aas,
                                bool permissive)
  {
    aas.peptide_.clear();
    String peptide(pep);
    peptide.trim();

    if (peptide.empty()) return;

    // remove optional n and c at start and end of string
    if (peptide[0] == 'n') 
    {
      peptide.erase(0,1);
    }

    if (peptide.empty()) return;

    if (peptide[peptide.size()-1] == 'c') 
    {
      peptide.erase(peptide.size()-1, 1);
    }

    if (peptide.empty()) return;
    
    // detect if this is the new dot notation containing dots for termini and
    // track if last char denoted a terminus
    bool dot_terminal(false), dot_notation(false);

    static ResidueDB* rdb = ResidueDB::getInstance();

    for (String::ConstIterator str_it = peptide.begin();
         str_it != peptide.end(); ++str_it)
    {
      // skip (optional) terminal delimiters
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
        aas.peptide_.push_back(r);
        continue;
      }

      // 2. modification: 
      //   determine specificity: 
      //     - at termini we first assume we are dealing with a N- or C-terminal modifications
      //       and fall back to (internal) modifications if there is none in our DB
      //     - otherwise we can be sure we are dealing with an internal modification
      ResidueModification::TermSpecificity specificity = ResidueModification::ANYWHERE;

      //   at the terminus we assume we are dealing with a N- or C-terminal modifications

      // make str_it point on '[' and set specificity if we are dealing with a terminus
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
        str_it = parseModRoundBrackets_(str_it, peptide, aas, dot_notation, dot_terminal);
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

    // deal with a single, unmodified X residue which is an indication for a
    // problem:
    // While PEPTIX[123]DE makes sense and represents an unknown mass of 123.0
    // Da, the sequence PEPTIXDE does not make sense as it is unclear what a
    // single, unknown residue should represent.
    const Residue* x = rdb->getResidue("X");
    if (aas.has(*x) )
    {
      std::cerr << "Error while parsing sequence " << pep << ": found an unknown AA without an estimated mass. Please use PEPTIX[123] syntax to indicate an unknown amino acid with a known mass." << std::endl;
      // throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error while parsing", String(pep));
    }
  }

  void AASequence::getAAFrequencies(Map<String, Size>& frequency_table) const
  {
    OPENMS_PRECONDITION( !this->has( *ResidueDB::getInstance()->getResidue("X") ), "Cannot get AA frequencies of sequence with unknown AA with unknown mass.");

    frequency_table.clear();

    for (std::vector<const Residue*>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
    {
      frequency_table[(*it)->getOneLetterCode()] += 1;
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

  void AASequence::setNTerminalModification(const String& modification)
  {
    if (modification == "")
    {
      n_term_mod_ = nullptr;
      return;
    }

    n_term_mod_ = &ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::N_TERM);
  }

  void AASequence::setCTerminalModification(const String& modification)
  {
    if (modification == "")
    {
      c_term_mod_ = nullptr;
      return;
    }
    c_term_mod_ = &ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::C_TERM);
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
