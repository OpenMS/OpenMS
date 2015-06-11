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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

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
    n_term_mod_(0),
    c_term_mod_(0)
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

  const Residue& AASequence::getResidue(SignedSize index) const
  {
    if (index >= 0 && Size(index) >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
    }
    if (index < 0)
    {
      throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
    }
    return *peptide_[index];
  }

  const Residue& AASequence::getResidue(Size index) const
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
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

    if (peptide_.size() >= 1)
    {
      // x/c-ion for a single residue is not defined
      if ((peptide_.size() == 1) && (type == Residue::XIon ||
       type == Residue::CIon))
      {
        LOG_ERROR << "AASequence::getFormula: Mass for ResidueType " << type << " not defined for sequences of length 1." << std::endl;
      }

      // Initialize with the missing/additional protons
      EmpiricalFormula ef; // = EmpiricalFormula("H") * charge; ??
      ef.setCharge(charge);

      // terminal modifications
      if (n_term_mod_ != 0 &&
        (type == Residue::Full || type == Residue::AIon ||
         type == Residue::BIon || type == Residue::CIon ||
         type == Residue::NTerminal)
        )
      {
        ef += n_term_mod_->getDiffFormula();
      }


      if (c_term_mod_ != 0 &&
        (type == Residue::Full || type == Residue::XIon ||
         type == Residue::YIon || type == Residue::ZIon ||
         type == Residue::CTerminal)
        )
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
        return ef + AASequence::getInternalToFull();

        case Residue::Internal:
        return ef;

        case Residue::NTerminal:
        return ef + AASequence::getInternalToNTerm();

        case Residue::CTerminal:
        return ef + AASequence::getInternalToCTerm();

        case Residue::AIon:
        return ef + AASequence::getInternalToAIon();

        case Residue::BIon:
        return ef + AASequence::getInternalToBIon();

        case Residue::CIon:
        return ef + AASequence::getInternalToCIon();

        case Residue::XIon:
        return ef + AASequence::getInternalToXIon();

        case Residue::YIon:
        return ef + AASequence::getInternalToYIon();

        case Residue::ZIon:
        return ef + AASequence::getInternalToZIon();

        default:
        LOG_ERROR << "AASequence::getFormula: unknown ResidueType" << std::endl;
      }

      return ef;
    }
    else
    {
      LOG_ERROR << "AASequence::getFormula: Formula for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
    }
  }

  double AASequence::getAverageWeight(Residue::ResidueType type, Int charge) const
  {
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
    
    if (peptide_.size() >= 1)
    {
      // x/c-ion for a single residue is not defined
      if ((peptide_.size() == 1) && (type == Residue::XIon ||
                                     type == Residue::CIon))
      {
        LOG_ERROR << "AASequence::getMonoWeight: Mass for ResidueType " << type << " not defined for sequences of length 1." << std::endl;
      }

      double mono_weight(Constants::PROTON_MASS_U * charge);

      // terminal modifications
      if (n_term_mod_ != 0 &&
          (type == Residue::Full || type == Residue::AIon ||
           type == Residue::BIon || type == Residue::CIon ||
           type == Residue::NTerminal))
      {
        mono_weight += n_term_mod_->getDiffMonoMass();
      }

      if (c_term_mod_ != 0 && 
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
        return mono_weight + AASequence::getInternalToFull().getMonoWeight();

        case Residue::Internal:
        return mono_weight;

        case Residue::NTerminal:
        return mono_weight + AASequence::getInternalToNTerm().getMonoWeight();

        case Residue::CTerminal:
        return mono_weight + AASequence::getInternalToCTerm().getMonoWeight();

        case Residue::AIon:
        return mono_weight + AASequence::getInternalToAIon().getMonoWeight();

        case Residue::BIon:
        return mono_weight + AASequence::getInternalToBIon().getMonoWeight();

        case Residue::CIon:
        return mono_weight + AASequence::getInternalToCIon().getMonoWeight();

        case Residue::XIon:
        return mono_weight + AASequence::getInternalToXIon().getMonoWeight();

        case Residue::YIon:
        return mono_weight + AASequence::getInternalToYIon().getMonoWeight();

        case Residue::ZIon:
        return mono_weight + AASequence::getInternalToZIon().getMonoWeight();

        default:
        LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
      }

      return mono_weight;
  }
  else
  {
    LOG_ERROR << "AASequence::getMonoWeight: Mass for ResidueType " << type << " not defined for sequences of length 0." << std::endl;
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

  const Residue& AASequence::operator[](SignedSize index) const
  {
    if (index < 0)
    {
      throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
    }
    else
    {
      if (Size(index) >= size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
      }
    }
    return *peptide_[Size(index)];
  }

  const Residue& AASequence::operator[](Size index) const
  {
    if (index >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
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
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "given residue");
    }
    AASequence seq = *this;
    seq += residue;
    return seq;
  }

  AASequence& AASequence::operator+=(const Residue* residue)
  {
    if (!ResidueDB::getInstance()->hasResidue(residue))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "given residue");
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
    }
    if (index + num > size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index + num, size());
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
    if (n_term_mod_ != 0 || c_term_mod_ != 0)
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

  bool AASequence::isModified(Size position) const
  {
    if (position >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, position, peptide_.size());
    }

    return peptide_[position]->isModified();
  }

  std::ostream& operator<<(std::ostream& os, const AASequence& peptide)
  {
    if (peptide.n_term_mod_ != 0)
    {
      os << "(" << peptide.n_term_mod_->getId() << ")";
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
          os << "[" << precisionWrapper(peptide.peptide_[i]->getMonoWeight()) << "]";
        }
        String id = ModificationsDB::getInstance()->getModification(
                      peptide.peptide_[i]->getOneLetterCode(),
                      peptide.peptide_[i]->getModification(),
                      ResidueModification::ANYWHERE).getId();
        if (id != "")
        {
          os << "(" << id << ")";
        }
        else
        {
          os << "([" << precisionWrapper(ModificationsDB::getInstance()->getModification(
                  peptide.peptide_[i]->getOneLetterCode(),
                  peptide.peptide_[i]->getModification(),
                  ResidueModification::ANYWHERE).getDiffMonoMass()) << "])";
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
            os << "[" << precisionWrapper(peptide.peptide_[i]->getMonoWeight()) << "]";
          }
        }
      }
    }

    if (peptide.c_term_mod_ != 0)
    {
      os << "(" << peptide.c_term_mod_->getId() << ")";
    }
    return os;
  }


  String::ConstIterator AASequence::parseModRoundBrackets_(
    const String::ConstIterator str_it, const String& str, AASequence& aas)
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
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, str,
          "Cannot convert string to peptide modification: missing ')'");
    }
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    if (aas.peptide_.empty()) // start of peptide -> N-terminal mod.
    {
      aas.n_term_mod_ = &(mod_db->getTerminalModification(
                            mod, ResidueModification::N_TERM));
      return mod_end;
    }
    if (std::distance(mod_end, str.end()) == 1) // end of peptide -> C-terminal mod.?
    {
      try
      {
        const ResidueModification* term_mod =
          &(mod_db->getTerminalModification(mod, ResidueModification::C_TERM));
        aas.c_term_mod_ = term_mod;
        return mod_end;
      }
      catch (Exception::ElementNotFound& /* e */)
      { // just do nothing, the mod is presumably a non-terminal one
      }
    }
    aas.peptide_.back() = ResidueDB::getInstance()->
      getModifiedResidue(aas.peptide_.back(), mod);
    // @TODO: if mod isn't found, "InvalidValue" is raised - catch it here?
    return mod_end;
  }

  String::ConstIterator AASequence::parseModSquareBrackets_(
    const String::ConstIterator str_it, const String& str, AASequence& aas)
  {
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
    String::ConstIterator mod_start = str_it;
    String::ConstIterator mod_end = ++mod_start;
    while ((mod_end != str.end()) && (*mod_end != ']')) ++mod_end;
    std::string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, str,
          "Cannot convert string to peptide modification: missing ']'");
    }

    double mass = String(mod).toDouble();
    size_t decimal_pos = mod.find('.');
    bool integer_mass = decimal_pos == std::string::npos;
    double tolerance = 0.5; // for integer mass values
    if (!integer_mass) // float mass values -> adapt tol. to decimal precision
    {
      size_t n_decimals = mod.size() - decimal_pos - 2;
      tolerance = std::pow(10.0, -int(n_decimals));
    }
    bool delta_mass = (mod[0] == '+') || (mod[0] == '-');
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    const Residue* residue = 0;
    if (!aas.peptide_.empty())
    {
      // internal modification (why not potentially C-terminal?):
      residue = aas.peptide_.back();
      if (delta_mass && (residue->getMonoWeight() <= 0.0)) // not allowed
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, str,
            "Using a mass difference to specify a modification on a residue of unknown mass is not supported in '" + \
            residue->getOneLetterCode() + "[" + mod +  "]'");
      }
      if (integer_mass) // use first modification that matches approximately
      {
        std::vector<String> res_mods;
        if (!delta_mass) // compute delta mass based on residue mass
        {
          // we expect that delta mass is relative to the full mass of the
          // residue, not its "internal" mass in the peptide (loss of H2O)!
          mass -= residue->getMonoWeight(Residue::Internal);
          delta_mass = true; // in case we need to create a new residue below
        }
        mod_db->getModificationsByDiffMonoMass(
          res_mods, residue->getOneLetterCode(), mass, tolerance);
        if (!res_mods.empty())
        {
          const ResidueModification* res_mod = &(mod_db->
                                                 getModification(res_mods[0]));
          aas.peptide_.back() = ResidueDB::getInstance()->
            getModifiedResidue(residue, res_mod->getId());
          return mod_end;
        }
      }
      else // float mass -> use best-matching modification
      {
        const ResidueModification* res_mod;
        if (delta_mass)
        {
          res_mod = mod_db->getBestModificationsByDiffMonoMass(
            residue->getOneLetterCode(), mass, tolerance);
        }
        else // absolute mass
        {
          res_mod = mod_db->getBestModificationsByMonoMass(
            residue->getOneLetterCode(), mass, tolerance);
        }
        if (res_mod)
        {
          aas.peptide_.back() = ResidueDB::getInstance()->
            getModifiedResidue(residue, res_mod->getId());
          return mod_end;
        }
      }
      LOG_WARN << "Warning: unknown modification '" + mod + "' of residue '" +
        residue->getOneLetterCode() + "' - adding it to the database" << std::endl;
    }
    // at beginning of peptide:
    else if (delta_mass) // N-terminal mod can only be specified by delta mass
    {
      std::vector<String> term_mods;
      mod_db->getTerminalModificationsByDiffMonoMass(
        term_mods, mass, tolerance, ResidueModification::N_TERM);
      if (!term_mods.empty())
      {
        aas.n_term_mod_ = &(mod_db->getTerminalModification(
                              term_mods[0], ResidueModification::N_TERM));
        return mod_end;
      }
      LOG_WARN << "Warning: unknown N-terminal modification '" + mod + "' - adding it to the database" << std::endl;
    }
    // create new modification:
    Residue new_res;
    new_res.setName(mod);
    if (residue && delta_mass)
    {
      new_res.setMonoWeight(mass + residue->getMonoWeight());
      new_res.setAverageWeight(mass + residue->getAverageWeight());
    }
    else
    { // mass value is for an internal residue, but methods expect full residue:
      new_res.setMonoWeight(mass + Residue::getInternalToFullMonoWeight());
      new_res.setAverageWeight(mass +
                               Residue::getInternalToFullAverageWeight());
    }
    ResidueDB::getInstance()->addResidue(new_res);
    aas.peptide_.back() = ResidueDB::getInstance()->getResidue(mod);
    return mod_end;
  }

  void AASequence::parseString_(const String& pep, AASequence& aas,
                                bool permissive)
  {
    aas.peptide_.clear();
    String peptide(pep);
    peptide.trim();
    if (peptide.empty()) return;

    static ResidueDB* rdb = ResidueDB::getInstance();
    for (String::ConstIterator str_it = peptide.begin();
         str_it != peptide.end(); ++str_it)
    {
      const Residue* r = rdb->getResidue(*str_it); // "isalpha" check not needed
      if (r)
      {
        aas.peptide_.push_back(r);
      }
      else if (*str_it == '(')
      {
        str_it = parseModRoundBrackets_(str_it, peptide, aas);
      }
      else if (*str_it == '[')
      {
        str_it = parseModSquareBrackets_(str_it, peptide, aas);
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
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide,
              "Cannot convert string to amino acid sequence: unexpected character '" + String(*str_it) + "'");
        }
      }
    }
  }

  void AASequence::getAAFrequencies(Map<String, Size>& frequency_table) const
  {
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
    }
    peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], modification);
  }

  void AASequence::setNTerminalModification(const String& modification)
  {
    if (modification == "")
    {
      n_term_mod_ = 0;
      return;
    }
    n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::N_TERM);
  }

  void AASequence::setCTerminalModification(const String& modification)
  {
    if (modification == "")
    {
      c_term_mod_ = 0;
      return;
    }
    c_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::C_TERM);
  }

  const String& AASequence::getNTerminalModification() const
  {
    static const String mod = "";
    if (n_term_mod_ == 0)
    {
      return mod;
    }
    return n_term_mod_->getId();
  }

  const String& AASequence::getCTerminalModification() const
  {
    static const String mod = "";
    if (c_term_mod_ == 0)
    {
      return mod;
    }
    return c_term_mod_->getId();
  }

  bool AASequence::hasNTerminalModification() const
  {
    return n_term_mod_ != 0;
  }

  bool AASequence::hasCTerminalModification() const
  {
    return c_term_mod_ != 0;
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
