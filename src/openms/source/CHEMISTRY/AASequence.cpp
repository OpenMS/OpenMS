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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <OpenMS/CONCEPT/PrecisionWrapper.h>

#include <algorithm>

#include <iostream>
#include <sstream>

using namespace std;

namespace OpenMS
{
  AASequence::AASequence() :
    n_term_mod_(0),
    c_term_mod_(0)
  {
  }

  AASequence::AASequence(const AASequence & rhs) :
    peptide_(rhs.peptide_),
    n_term_mod_(rhs.n_term_mod_),
    c_term_mod_(rhs.c_term_mod_)
  {
  }

  AASequence::~AASequence()
  {
  }

  AASequence & AASequence::operator=(const AASequence & rhs)
  {
    if (this != &rhs)
    {
      peptide_ = rhs.peptide_;
      n_term_mod_ = rhs.n_term_mod_;
      c_term_mod_ = rhs.c_term_mod_;
    }
    return *this;
  }

  const Residue & AASequence::getResidue(SignedSize index) const
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

  const Residue & AASequence::getResidue(Size index) const
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
    }
    return *peptide_[index];
  }

  String AASequence::toString() const
  {
    stringstream ss;
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

  bool AASequence::operator<(const AASequence & rhs) const
  {
    return toString() < rhs.toString();
  }

  EmpiricalFormula AASequence::getFormula(Residue::ResidueType type, Int charge) const
  {
    EmpiricalFormula ef;
    ef.setCharge(charge);
    static EmpiricalFormula H("H");
    static EmpiricalFormula OH("OH");
    static EmpiricalFormula NH("NH");

    // terminal modifications
    if (n_term_mod_ != 0 &&
        (type == Residue::Full || type == Residue::AIon || type == Residue::BIon || type == Residue::CIon || type == Residue::NTerminal)
        )
    {
      ef += n_term_mod_->getDiffFormula();
    }


    if (c_term_mod_ != 0 &&
        (type == Residue::Full || type == Residue::XIon || type == Residue::YIon || type == Residue::ZIon || type == Residue::CTerminal)
        )
    {
      ef += c_term_mod_->getDiffFormula();
    }

    if (peptide_.size() > 0)
    {
      if (peptide_.size() == 1)
      {
        ef += peptide_[0]->getFormula(type);
      } else
      {
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
          return ef /* + add_protons*/;

        case Residue::NTerminal:
          return ef + Residue::getInternalToFull() - Residue::getNTerminalToFull();

        case Residue::CTerminal:
          return ef + Residue::getInternalToFull() - Residue::getCTerminalToFull();

        case Residue::BIon:
          return ef + Residue::getInternalToFull() - Residue::getBIonToFull() - H;

        case Residue::AIon:
          return ef + Residue::getInternalToFull() - Residue::getAIonToFull() - H;

        case Residue::CIon:
          return ef + Residue::getInternalToFull() - OH + NH;

        case Residue::XIon:
          return ef + Residue::getInternalToFull() + Residue::getXIonToFull();

        case Residue::YIon:
          return ef + Residue::getInternalToFull() + Residue::getYIonToFull();

        case Residue::ZIon:
          return ef + Residue::getInternalToFull() - Residue::getZIonToFull();

        default:
          cerr << "AASequence::getFormula: unknown ResidueType" << endl;
        }
      }
    }

    return ef;
  }

  double AASequence::getAverageWeight(Residue::ResidueType type, Int charge) const
  {
    // check whether tags are present
    double tag_offset(0);
    for (ConstIterator it = this->begin(); it != this->end(); ++it)
    {
      if (it->getOneLetterCode() == "")
      {
        tag_offset += it->getMonoWeight();
      }
    }
    return tag_offset + getFormula(type, charge).getAverageWeight();
  }

  double AASequence::getMonoWeight(Residue::ResidueType type, Int charge) const
  {
    // check whether tags are present
    double tag_offset(0);
    for (ConstIterator it = this->begin(); it != this->end(); ++it)
    {
      if (it->getOneLetterCode() == "")
      {
        tag_offset += it->getMonoWeight();
      }
    }

    return tag_offset + getFormula(type, charge).getMonoWeight();
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

  const Residue & AASequence::operator[](SignedSize index) const
  {
    if (index < 0)
    {
      throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, 0);
    } else
    {
      if (Size(index) >= size())
      {
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
      }
    }
    return *peptide_[Size(index)];
  }

  const Residue & AASequence::operator[](Size index) const
  {
    if (index >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, size());
    }
    return *peptide_[index];
  }

  AASequence & AASequence::operator+=(const AASequence & sequence)
  {
    for (Size i = 0; i != sequence.peptide_.size(); ++i)
    {
      peptide_.push_back(sequence.peptide_[i]);
    }
    return *this;
  }

  AASequence AASequence::operator+(const AASequence & sequence) const
  {
    AASequence seq;
    seq.peptide_ = peptide_;
    for (Size i = 0; i != sequence.peptide_.size(); ++i)
    {
      seq.peptide_.push_back(sequence.peptide_[i]);
    }
    return seq;
  }

  AASequence AASequence::operator+(const Residue * residue) const
  {
    if (!ResidueDB::getInstance()->hasResidue(residue))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "given residue");
    }
    AASequence seq = *this;
    seq += residue;
    return seq;
  }

  AASequence & AASequence::operator+=(const Residue * residue)
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
    for (Size i = 0; i < index; ++i)
    {
      seq.peptide_.push_back(peptide_[i]);
    }
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
    for (Size i = size() - index; i != size(); ++i)
    {
      seq.peptide_.push_back(peptide_[i]);
    }
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
    for (Size i = index; i != index + num; ++i)
    {
      seq.peptide_.push_back(peptide_[i]);
    }
    return seq;
  }

  bool AASequence::has(const Residue & residue) const
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

  bool AASequence::hasSubsequence(const AASequence & sequence) const
  {
    if (sequence.empty())
    {
      return true;
    } else
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
              } else
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

  bool AASequence::hasPrefix(const AASequence & sequence) const
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

  bool AASequence::hasSuffix(const AASequence & sequence) const
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

  bool AASequence::operator==(const AASequence & peptide) const
  {
    if (size() != peptide.size())
    {
      return false;
    }

    for (Size i = 0; i != size(); ++i)
    {
      if (peptide_[i] != peptide.peptide_[i])
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

  bool AASequence::operator!=(const AASequence & peptide) const
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

    for (vector<const Residue *>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
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
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide_.size(), position);
    }

    return peptide_[position]->isModified();
  }

  ostream & operator<<(ostream & os, const AASequence & peptide)
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
        } else
        {
          os << "[" << precisionWrapper(peptide.peptide_[i]->getMonoWeight()) << "]";
        }
        String id = ModificationsDB::getInstance()->getModification(peptide.peptide_[i]->getOneLetterCode(), peptide.peptide_[i]->getModification(), ResidueModification::ANYWHERE).getId();
        if (id != "")
        {
          os << "(" << id << ")";
        } else
        {
          os << "([" << precisionWrapper(ModificationsDB::getInstance()->getModification(peptide.peptide_[i]->getOneLetterCode(), peptide.peptide_[i]->getModification(), ResidueModification::ANYWHERE).getDiffMonoMass()) << "])";
        }
      } else
      {
        if (peptide.peptide_[i]->getOneLetterCode() != "")
        {
          os << peptide.peptide_[i]->getOneLetterCode();
        } else
        {
          if (peptide.peptide_[i]->getShortName() != "")
          {
            os << peptide.peptide_[i]->getShortName();
          } else
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

  void AASequence::parseString_(AASequence& aas, const String & pep)
  {
    aas.peptide_.clear();
    String peptide(pep);
    peptide.trim();

    if (peptide.empty())
    {
      return;
    }

    // split the peptide in its residues
    vector<String> split;
    Size pos(0);
    bool mod_open(false);

    if (peptide[0] == '(')
    {
      mod_open = true;
    }
    Size num_brackets(0);
    const Size& size_peptide = peptide.size();
    for (Size i = 1; i < size_peptide; ++i)
    {
      if (isalpha(peptide[i]) && isupper(peptide[i]) && !mod_open)
      {
        split.push_back(peptide.substr(pos, i - pos));
        pos = i;
      }
      switch (peptide[i])
      {
        case '(':
          if (mod_open)
          {
            ++num_brackets;
            continue;
          }
          mod_open = true;
          break;
        case ')':
          if (num_brackets != 0)
          {
            --num_brackets;
            continue;
          }
          mod_open = false;
          break;
        case '[':
          if (mod_open)
          {
            ++num_brackets;
            continue;
          }
          mod_open = true;
          break;
        case ']':
          if (num_brackets != 0)
          {
            --num_brackets;
            continue;
          }
          mod_open = false;
          break;
        default:
          break;
      }
    }

    // push_back last residue
    split.push_back(peptide.substr(pos, peptide.size() - pos));

    if (!split.empty() && !split[0].empty() && split[0][0] == '(')
    {
      String mod = split[0];
      mod.trim();
      mod.erase(mod.begin());
      mod.erase(mod.end() - 1);
      aas.n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(mod, ResidueModification::N_TERM);

      split.erase(split.begin());
    }

    if (split.empty())
    {
      return;
    }

    // test the last split if there is a C-terminal modification
    const String c_term = *(split.end() - 1);
    Size c_term_mods = count(c_term.begin(), c_term.end(), '(');
    if (c_term_mods > 0)
    {
      // now we have found a potential C-term modification
      String mod;

      // we need to parse here since K(Label:18O(2)) .. will not be found
      // correctly by prefix/suffix search
      Size brackets = 0;
      // we start at (end - 1) to skip trailing ')'
      for (String::ConstReverseIterator rIt = (c_term.rbegin() + 1); rIt != c_term.rend(); ++rIt)
      {
        if (*rIt == '(' && brackets == 0)
        {
          break; // we reached beginning of mod
        } else if (*rIt == ')')
        {
          ++brackets;
        } else if (*rIt == '(')
        {
          --brackets;
        }

        mod = *rIt + mod;
      }

      try
      {
        const ResidueModification * potential_mod = &ModificationsDB::getInstance()->getTerminalModification(mod, ResidueModification::C_TERM);
        aas.c_term_mod_ = potential_mod;
        split[split.size() - 1] = c_term.substr(0, c_term.size() - mod.size() - 2);
      }
      catch (Exception::ElementNotFound & /* e */)
      {
        // just do nothing the mod is presumably a non-terminal one
      }
    }

    // parse the residues
    for (Size i = 0; i != split.size(); ++i)
    {
      const String& res = split[i];
      String name, mod, tag;
      for (Size j = 0; j != res.size(); ++j)
      {
        if (isalpha(res[j]))
        {
          name += res[j];
        } else
        {
          if (res[j] == '(')
          {
            if (res[res.size() - 1] != ')')
            {
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Missing ')'!");
            }

            for (Size k = j + 1; k < res.size() - 1; ++k)             // skip last ')'
            {
              mod += res[k];
            }
            break;
          } else
          {
            if (res[j] == '[')
            {
              for (Size k = j + 1; res[k] != ']'; ++k)
              {
                if (k == res.size())
                {
                  throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Missing ']'!");
                }
                tag += res[k];
              }
              break;
            } else
            {
              throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Residue '" + res + "' unknown at position " + String(j) + ", residue # " + String(i) + " !");
            }
          }
        }
      }

      // Retrieve the underlying residue
      const Residue * res_ptr = ResidueDB::getInstance()->getResidue(name);

      if (res_ptr == 0 && tag.empty() && mod.empty())
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Cannot parse residue with name: '" + name + "'!");
        return;
      }

      if (!mod.empty() && i > 0)
      {
        if (res_ptr == 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Cannot parse residue with name: '" + name + "'!");
        }
        aas.peptide_.push_back(ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod));
      } else if (!mod.empty() && i == 0)
      {
        if (res_ptr == 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Cannot parse residue with name: '" + name + "'!");
        }
        std::set<const ResidueModification *> mod_candidates;
        ModificationsDB::getInstance()->searchModifications(mod_candidates, res_ptr->getOneLetterCode(), mod, ResidueModification::ANYWHERE);
        if (!mod_candidates.empty())
        {
          aas.peptide_.push_back(ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod));
        } else
        {
          ModificationsDB::getInstance()->searchTerminalModifications(mod_candidates, mod, ResidueModification::N_TERM);
          if (!mod_candidates.empty())
          {
            aas.n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(mod, ResidueModification::N_TERM);
            aas.peptide_.push_back(res_ptr);
          }
        }
      } else if (!tag.empty())
      {
        bool is_NTerm = (res_ptr == 0) && name.empty() ? true : false;

        if (tag.hasPrefix("+") || tag.hasPrefix("-"))
        {
          // delta mass
          double delta_mass = tag.toDouble();
          const Residue* result = NULL;

          if (tag.hasSubstring("."))
          {
            // signed float tag [+123.456] -> look for an exact match
            if (is_NTerm)
            {
              vector<String> mods;
              ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, delta_mass, 1.0, ResidueModification::N_TERM);
              if (!mods.empty())
              {
                aas.n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(mods[0], ResidueModification::N_TERM);
                continue;
              }
            } else
            {
              if (res_ptr == 0)
              {
                throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Cannot parse residue with name: '" + name + "'!");
              }
              const ResidueModification * mod = ModificationsDB::getInstance()->getBestModificationsByDiffMonoMass(res_ptr->getOneLetterCode(), delta_mass, 1.0);
              if (mod != NULL)
              {
                result = ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod->getId());
              }
            }
          } else
          {
            // signed integer tag [+123] -> look for the first match (usually this is what is intended)
            if (is_NTerm)
            {
              vector<String> mods;
              ModificationsDB::getInstance()->getTerminalModificationsByDiffMonoMass(mods, delta_mass, 0.5, ResidueModification::N_TERM);
              if (!mods.empty())
              {
                aas.n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(mods[0], ResidueModification::N_TERM);
                continue;
              }
            } else
            {
              std::vector< String > mods;
              if (res_ptr == 0)
              {
                throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Cannot parse residue with name: '" + name + "'!");
              }
              ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, res_ptr->getOneLetterCode(), delta_mass, 0.5);
              if (!mods.empty())
              {
                const ResidueModification * mod = &ModificationsDB::getInstance()->getModification(mods[0]);
                result = ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod->getId());
              }
            }
          }

          // using an amino acid with zero weight and a differential modification on that cannot lead to a valid result!
          if (res_ptr && res_ptr->getMonoWeight() <= 0.0)
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Having a difference modification on an unspecified residue (" + name + "[" + tag +  "]) will probably not produce a correct mass.");
          }

          if (result == NULL)
          {
            std::cout <<  "Warning: unknown modification " << tag << " on residue "
              << name << ": will add to the database." << std::endl;
            Residue res(tag, String(""), String(""), EmpiricalFormula(""));
            if (is_NTerm)
            {
              res.setMonoWeight(delta_mass);
              res.setAverageWeight(delta_mass);
            } else
            {
              res.setMonoWeight(delta_mass + res_ptr->getMonoWeight());
              res.setAverageWeight(delta_mass + res_ptr->getAverageWeight());
            }
            ResidueDB::getInstance()->addResidue(res);
            result = ResidueDB::getInstance()->getResidue(tag);
          }
          aas.peptide_.push_back(result);
        } else  // absolute (unsigned) mass tags []
        {
          double mass = tag.toDouble();
          const Residue* result = NULL;

          if (tag.hasSubstring("."))
          {
            // we have a float, look for an exact match
            const ResidueModification * mod = ModificationsDB::getInstance()->getBestModificationsByMonoMass(res_ptr->getOneLetterCode(), mass, 1.0);
            if (mod)
              result = ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod->getId());
          } else
          {
            // we have an integer, look for the first match (usually this is what is intended)
            std::vector< String > mods;
            double res_deltamass = mass - (res_ptr->getMonoWeight() - res_ptr->getInternalToFullMonoWeight());
            ModificationsDB::getInstance()->getModificationsByDiffMonoMass(mods, res_ptr->getOneLetterCode(), res_deltamass, 0.5);
            if (!mods.empty())
            {
              const ResidueModification * mod = &ModificationsDB::getInstance()->getModification(mods[0]);
              result = ResidueDB::getInstance()->getModifiedResidue(res_ptr, mod->getId());
            }
          }

          if (result == NULL || res_ptr->getMonoWeight() <= 0.0)
          {
            // using an amino acid with zero weight should lead to an accurate mass representation of the AA (and not try to guess something)... !
            std::cout <<  "Warning: unknown modification " << tag << " on residue " << name << ": will add to the database." << std::endl;
            Residue res(tag, String(""), String(""), EmpiricalFormula(""));
            res.setMonoWeight(mass);
            res.setAverageWeight(mass);
            ResidueDB::getInstance()->addResidue(res);
            result = ResidueDB::getInstance()->getResidue(tag);
          }
          aas.peptide_.push_back(result);
        }
      } else
      {
        // mod and tag are both empty
        aas.peptide_.push_back(res_ptr);
        if (aas.peptide_.size() < i)
        {
          throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, peptide, "Cannot convert string into AASequence. Mod and tag are both empty.");
          return;
        }
      }
    }
  }

  void AASequence::getAAFrequencies(Map<String, Size> & frequency_table) const
  {
    frequency_table.clear();

    for (vector<const Residue *>::const_iterator it = peptide_.begin(); it != peptide_.end(); ++it)
    {
      frequency_table[(*it)->getOneLetterCode()] += 1;
    }
  }

  void AASequence::setModification(Size index, const String & modification)
  {
    if (index >= peptide_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, peptide_.size());
    }
    peptide_[index] = ResidueDB::getInstance()->getModifiedResidue(peptide_[index], modification);
  }

  void AASequence::setNTerminalModification(const String & modification)
  {
    if (modification == "")
    {
      n_term_mod_ = 0;
      return;
    }
    n_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::N_TERM);
  }

  void AASequence::setCTerminalModification(const String & modification)
  {
    if (modification == "")
    {
      c_term_mod_ = 0;
      return;
    }
    c_term_mod_ = &ModificationsDB::getInstance()->getTerminalModification(modification, ResidueModification::C_TERM);
  }

  const String & AASequence::getNTerminalModification() const
  {
    static const String mod = "";
    if (n_term_mod_ == 0)
    {
      return mod;
    }
    return n_term_mod_->getId();
  }

  const String & AASequence::getCTerminalModification() const
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

  AASequence AASequence::fromString(const String & s)
  {
    AASequence aas;
    parseString_(aas, s);
    return aas;
  }

  AASequence AASequence::fromString(const char * s)
  {
    AASequence aas;
    parseString_(aas, String(s));
    return aas;
  }
}
