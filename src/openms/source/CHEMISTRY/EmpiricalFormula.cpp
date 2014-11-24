// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  EmpiricalFormula::EmpiricalFormula() :
    charge_(0)
  {
  }

  EmpiricalFormula::EmpiricalFormula(const EmpiricalFormula& formula) :
    formula_(formula.formula_),
    charge_(formula.charge_)
  {
  }

  EmpiricalFormula::EmpiricalFormula(const String& formula)
  {
    charge_ = parseFormula_(formula_, formula);
  }

  EmpiricalFormula::EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge)
  {
    formula_[element] = number;
    charge_ = charge;
  }

  EmpiricalFormula::~EmpiricalFormula()
  {
  }

  double EmpiricalFormula::getMonoWeight() const
  {
    double weight(0);
    if (charge_ > 0)
    {
      weight += Constants::PROTON_MASS_U * charge_;
    }
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      weight += it->first->getMonoWeight() * (double)it->second;
    }
    return weight;
  }

  double EmpiricalFormula::getAverageWeight() const
  {
    double weight(0);
    if (charge_ > 0)
    {
      weight += Constants::PROTON_MASS_U * charge_;
    }
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      weight += it->first->getAverageWeight() * (double)it->second;
    }
    return weight;
  }

  IsotopeDistribution EmpiricalFormula::getIsotopeDistribution(UInt max_depth) const
  {
    IsotopeDistribution result(max_depth);
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      IsotopeDistribution tmp = it->first->getIsotopeDistribution();
      tmp.setMaxIsotope(max_depth);
      result += tmp * it->second;
    }
    result.renormalize();
    return result;
  }

  SignedSize EmpiricalFormula::getNumberOf(const Element* element) const
  {
    MapType_::const_iterator it  = formula_.find(element);
    if (it != formula_.end())
    {
      return it->second;
    }
    return 0;
  }

  SignedSize EmpiricalFormula::getNumberOfAtoms() const
  {
    SignedSize num_atoms(0);
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      num_atoms += it->second;
    }
    return num_atoms;
  }

  void EmpiricalFormula::setCharge(SignedSize charge)
  {
    charge_ = charge;
  }

  SignedSize EmpiricalFormula::getCharge() const
  {
    return charge_;
  }

  String EmpiricalFormula::toString() const
  {
    String formula;
    std::map<String, SignedSize> new_formula;

    for (MapType_::const_iterator it = formula_.begin(); it != formula_.end(); ++it)
    {
      new_formula[it->first->getSymbol()] = it->second;
    }

    for (std::map<String, SignedSize>::const_iterator it = new_formula.begin(); it != new_formula.end(); ++it)
    {
      formula += it->first + String(it->second);
    }
    return formula;
  }

  EmpiricalFormula& EmpiricalFormula::operator=(const EmpiricalFormula& formula)
  {
    if (this != &formula)
    {
      formula_ = formula.formula_;
      charge_ = formula.charge_;
    }
    return *this;
  }

  EmpiricalFormula EmpiricalFormula::operator*(const SignedSize& times) const
  {
    EmpiricalFormula ef(*this);
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      ef.formula_[it->first] *= times;
    }
    ef.charge_ *= times;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula EmpiricalFormula::operator+(const EmpiricalFormula& formula) const
  {
    EmpiricalFormula ef;
    ef.formula_ = formula.formula_;
    MapType_::const_iterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      MapType_::iterator ef_it  = ef.formula_.find(it->first);
      if (ef_it != ef.formula_.end())
      {
        ef_it->second += it->second;
      }
      else
      {
        ef.formula_.insert(*it);
      }
    }
    ef.charge_ = charge_ + formula.charge_;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula& EmpiricalFormula::operator+=(const EmpiricalFormula& formula)
  {
    MapType_::const_iterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      MapType_::iterator f_it  = formula_.find(it->first);
      if (f_it != formula_.end())
      {
        f_it->second += it->second;
      }
      else
      {
        formula_.insert(*it);
      }
    }
    charge_ += formula.charge_;
    removeZeroedElements_();
    return *this;
  }

  EmpiricalFormula EmpiricalFormula::operator-(const EmpiricalFormula& formula) const
  {
    EmpiricalFormula ef(*this);
    MapType_::const_iterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      const Element* e = it->first;
      SignedSize num = it->second;
      MapType_::iterator ef_it  = ef.formula_.find(e);
      if (ef_it != ef.formula_.end())
      {
        ef_it->second -= num;
      }
      else
      {
        ef.formula_[e] = -num;
      }
    }

    ef.charge_ = charge_ - formula.charge_;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula& EmpiricalFormula::operator-=(const EmpiricalFormula& formula)
  {
    MapType_::const_iterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      MapType_::iterator f_it  = formula_.find(it->first);
      if (f_it != formula_.end())
      {
        f_it->second -= it->second;
      }
      else
      {
        formula_[it->first] = -it->second;
      }
    }
    charge_ -= formula.charge_;
    removeZeroedElements_();
    return *this;
  }

  bool EmpiricalFormula::isCharged() const
  {
    return charge_ != 0;
  }

  bool EmpiricalFormula::isEmpty() const
  {
    return formula_.empty();
  }

  bool EmpiricalFormula::hasElement(const Element* element) const
  {
    return formula_.find(element) != formula_.end();
  }

  bool EmpiricalFormula::contains(const EmpiricalFormula& ef)
  {
    for (EmpiricalFormula::const_iterator it = ef.begin(); it != ef.end(); ++it)
    {
      if (this->getNumberOf(it->first) < it->second)
      {
        return false;
      }
    }
    return true;
  }

  bool EmpiricalFormula::operator==(const EmpiricalFormula& formula) const
  {
    return formula_ == formula.formula_ && charge_ == formula.charge_;
  }

  bool EmpiricalFormula::operator!=(const EmpiricalFormula& formula) const
  {
    return formula_ != formula.formula_ || charge_ != formula.charge_;
  }

  ostream& operator<<(ostream& os, const EmpiricalFormula& formula)
  {
    std::map<String, SignedSize> new_formula;
    for (Map<const Element*, SignedSize>::const_iterator it = formula.formula_.begin(); it != formula.formula_.end(); ++it)
    {
      new_formula[it->first->getSymbol()] = it->second;
    }

    for (std::map<String, SignedSize>::const_iterator it = new_formula.begin(); it != new_formula.end(); ++it)
    {
      os << it->first;
      if (it->second > 1)
      {
        os << it->second;
      }
    }
    if (formula.charge_ == 0)
    {
      return os;
    }

    if (formula.charge_ > 0)
    {
      if (formula.charge_ == 1)
      {
        os << "+";
      }
      else
      {
        os << "+" << formula.charge_;
      }
    }
    else
    {
      if (formula.charge_ == -1)
      {
        os << "-";
      }
      else
      {
        os << "-" << formula.charge_;
      }
    }
    return os;
  }

  SignedSize EmpiricalFormula::parseFormula_(std::map<const Element*, SignedSize>& ef, const String& input_formula) const
  {
    SignedSize charge = 0;
    String formula(input_formula);

    // we start with the charge part, read until the begin of the formula or a element symbol occurs
    String suffix;
    for (SignedSize reverse_i(formula.size() - 1); reverse_i >= 0; --reverse_i)
    {
      if (!isalpha(formula[reverse_i]))
      {
        suffix = formula[reverse_i] + suffix;
      }
      else
      {
        break;
      }
    }

    // determine charge
    if (!suffix.empty())
    {
      Size i = 1;
      for (; i < suffix.size(); ++i)
      {
        if (!isdigit(suffix[i]))
        {
          break;
        }
      }
      if (i != suffix.size())
      {
        // we found the charge part
        String charge_str;
        for (Size j = i + 1; j < suffix.size(); ++j)
        {
          charge_str += suffix[j];
        }

        SignedSize tmp_charge = 1;
        if (!charge_str.empty())
        {
          tmp_charge = charge_str.toInt();
        }
        if (suffix[i] == '-')
        {
          charge = -1 * tmp_charge;
        }
        else
        {
          if (suffix[i] == '+')
          {
            charge = tmp_charge;
          }
          else
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, formula, "Cannot parse charge part of formula!");
          }
        }

        // now remove the charge part from the formula
        formula.resize(formula.size() - charge_str.size() - 1);
      }
    }

    if (suffix.size() == 1 && suffix[0] == '+')
    {
      charge = 1;
      formula.resize(formula.size() - 1);
    }
    else if (suffix.size() == formula.size())
    {
      if (suffix.size() > 1)
      {
        if (suffix[0] == '-' || suffix[0] == '+')
        {
          charge = suffix.toInt();
          return charge;
        }
      }
      else
      {
        if (suffix == "-")
        {
          charge = -1;
          return charge;
        }
      }
    }

    // split the formula
    vector<String> splitter;
    if (formula.size() > 0)
    {
      if (!isdigit(formula[0]) || formula[0] == '(')
      {
        bool is_isotope(false), is_symbol(false);
        String split;
        for (Size i = 0; i < formula.size(); ++i)
        {
          if ((isupper(formula[i]) && (!is_isotope || is_symbol))
             || formula[i] == '(')
          {
            if (split != "")
            {
              splitter.push_back(split);
              is_isotope = false;
              is_symbol = false;
            }
            split = String(1, formula[i]);
          }
          else
          {
            split += String(1, formula[i]);
          }
          if (formula[i] == '(')
          {
            is_isotope = true;
          }
          if (isupper(formula[i]))
          {
            is_symbol = true;
          }
        }
        splitter.push_back(split);
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, formula, "This formula does not begin with an element!");
      }
    }

    // add up the elements
    for (Size i = 0; i != splitter.size(); ++i)
    {
      String split = splitter[i];
      String number;
      String symbol;
      bool had_symbol(false);
      for (SignedSize j = split.size() - 1; j >= 0; --j)
      {
        if (!had_symbol && (isdigit(split[j]) || split[j] == '-'))
        {
          number = split[j] + number;
        }
        else
        {
          symbol = split[j] + symbol;
          had_symbol = true;
        }
      }

      SignedSize num(1);
      if (number != "")
      {
        num = number.toInt();
      }

      const ElementDB* db = ElementDB::getInstance();
      if (db->hasElement(symbol))
      {
        if (num != 0)
        {
          const Element* e = db->getElement(symbol);
          std::map<const Element*, SignedSize>::iterator it = ef.find(e);
          if (it != ef.end())
          {
            it->second += num;
          }
          else
          {
            ef.insert(std::make_pair(e, num));
          }
        }
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown element '" + split + "'", "'" + symbol + "' found. Please use only valid element identifiers or modify share/OpenMS/CHEMISTRY/Elements.xml!");
      }
    }

    // remove elements with 0 counts
    std::map<const Element*, SignedSize>::iterator it = ef.begin();
    while (it != ef.end())
    {
      if (it->second == 0)
      {
        ef.erase(it++);   // Note: post increment needed! Otherwise iterator is invalidated
      }
      else
      {
        ++it;
      }
    }

    return charge;
  }

  void EmpiricalFormula::removeZeroedElements_()
  {
    MapType_::iterator it = formula_.begin();
    while (it != formula_.end())
    {
      if (it->second == 0)
      {
        formula_.erase(it++);   // Note: post increment needed! Otherwise iterator is invalidated
      }
      else
      {
        ++it;
      }
    }
  }

}
