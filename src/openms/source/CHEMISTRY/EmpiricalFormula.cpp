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

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace std;

namespace OpenMS
{
  EmpiricalFormula::EmpiricalFormula() :
    charge_(0),
    element_db_(ElementDB::getInstance())
  {
  }

  EmpiricalFormula::EmpiricalFormula(const EmpiricalFormula & formula) :
    formula_(formula.formula_),
    charge_(formula.charge_),
    element_db_(formula.element_db_)
  {
  }

  EmpiricalFormula::EmpiricalFormula(const String & formula) :
    element_db_(ElementDB::getInstance())
  {
    charge_ = parseFormula_(formula_, formula);
  }

  EmpiricalFormula::EmpiricalFormula(SignedSize number, const Element * element, SignedSize charge)
  {
    formula_[element] = number;
    charge_ = charge;
  }

  EmpiricalFormula::~EmpiricalFormula()
  {
  }

  DoubleReal EmpiricalFormula::getMonoWeight() const
  {
    DoubleReal weight(0);
    if (charge_ > 0)
    {
      weight += Constants::PROTON_MASS_U * charge_;
    }
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      weight += it->first->getMonoWeight() * (DoubleReal)it->second;
    }
    return weight;
  }

  DoubleReal EmpiricalFormula::getAverageWeight() const
  {
    DoubleReal weight(0);
    if (charge_ > 0)
    {
      weight += Constants::PROTON_MASS_U * charge_;
    }
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      weight += it->first->getAverageWeight() * (DoubleReal)it->second;
    }
    return weight;
  }

  IsotopeDistribution EmpiricalFormula::getIsotopeDistribution(UInt max_depth) const
  {
    IsotopeDistribution result(max_depth);
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      IsotopeDistribution tmp = it->first->getIsotopeDistribution();
      tmp.setMaxIsotope(max_depth);
      result += tmp * it->second;
    }
    result.renormalize();
    return result;
  }

  const Element * EmpiricalFormula::getElement(const String & name) const
  {
    return element_db_->getElement(name);
  }

  const Element * EmpiricalFormula::getElement(UInt atomic_number) const
  {
    return element_db_->getElement(atomic_number);
  }

  const ElementDB * EmpiricalFormula::getElementDB() const
  {
    return element_db_;
  }

  Size EmpiricalFormula::getNumberOf(UInt atomic_number) const
  {
    if (element_db_->hasElement(atomic_number))
    {
      if (formula_.has(element_db_->getElement(atomic_number)))
      {
        return formula_[element_db_->getElement(atomic_number)];
      }
    }
    return 0;
  }

  Size EmpiricalFormula::getNumberOf(const String & name) const
  {
    if (element_db_->hasElement(name))
    {
      if (formula_.has(element_db_->getElement(name)))
      {
        return formula_[element_db_->getElement(name)];
      }
    }
    return 0;
  }

  Size EmpiricalFormula::getNumberOf(const Element * element) const
  {
    if (formula_.has(element))
    {
      return formula_[element];
    }
    return 0;
  }

  Size EmpiricalFormula::getNumberOfAtoms() const
  {
    Size num_atoms(0);
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
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
    Map<String, SignedSize> new_formula;

    for (Map<const Element *, SignedSize>::ConstIterator it = formula_.begin(); it != formula_.end(); ++it)
    {
      new_formula[it->first->getSymbol()] = it->second;
    }

    for (Map<String, SignedSize>::ConstIterator it = new_formula.begin(); it != new_formula.end(); ++it)
    {
      formula += it->first + String(it->second);
    }
    return formula;
  }

  EmpiricalFormula & EmpiricalFormula::operator=(const EmpiricalFormula & formula)
  {
    if (this != &formula)
    {
      formula_ = formula.formula_;
      charge_ = formula.charge_;
    }
    return *this;
  }

  EmpiricalFormula & EmpiricalFormula::operator=(const String & formula)
  {
    charge_ = 0;
    formula_.clear();
    charge_ = parseFormula_(formula_, formula);
    return *this;
  }

  EmpiricalFormula EmpiricalFormula::operator*(const SignedSize & times) const
  {
    EmpiricalFormula ef(*this);
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      ef.formula_[it->first] *= times;
    }
    ef.charge_ *= times;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula EmpiricalFormula::operator+(const EmpiricalFormula & formula) const
  {
    EmpiricalFormula ef;
    ef.formula_ = formula.formula_;
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      if (ef.formula_.has(it->first))
      {
        ef.formula_[it->first] += it->second;
      }
      else
      {
        ef.formula_[it->first] = it->second;
      }
    }
    ef.charge_ = charge_ + formula.charge_;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula EmpiricalFormula::operator+(const String & formula) const
  {
    EmpiricalFormula ef;
    SignedSize charge = parseFormula_(ef.formula_, formula);
    Map<const Element *, SignedSize>::ConstIterator it = formula_.begin();
    for (; it != formula_.end(); ++it)
    {
      if (ef.formula_.has(it->first))
      {
        ef.formula_[it->first] += it->second;
      }
      else
      {
        ef.formula_[it->first] = it->second;
      }
    }
    ef.charge_ = charge_ + charge;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula & EmpiricalFormula::operator+=(const EmpiricalFormula & formula)
  {
    Map<const Element *, SignedSize>::ConstIterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      if (formula_.has(it->first))
      {
        formula_[it->first] += it->second;
      }
      else
      {
        formula_[it->first] = it->second;
      }
    }
    charge_ += formula.charge_;
    removeZeroedElements_();
    return *this;
  }

  EmpiricalFormula & EmpiricalFormula::operator+=(const String & formula)
  {
    Map<const Element *, SignedSize> str_formula;
    SignedSize charge = parseFormula_(str_formula, formula);
    charge_ += charge;
    Map<const Element *, SignedSize>::ConstIterator it;
    for (it = str_formula.begin(); it != str_formula.end(); ++it)
    {
      if (formula_.has(it->first))
      {
        formula_[it->first] += it->second;
      }
      else
      {
        formula_[it->first] = it->second;
      }
    }
    removeZeroedElements_();
    return *this;
  }

  EmpiricalFormula EmpiricalFormula::operator-(const EmpiricalFormula & formula) const
  {
    EmpiricalFormula ef(*this);
    Map<const Element *, SignedSize>::ConstIterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      const Element * e = it->first;
      SignedSize num = it->second;
      if (formula_.has(e))
      {
        ef.formula_[e] -= num;
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

  EmpiricalFormula EmpiricalFormula::operator-(const String & formula) const
  {
    EmpiricalFormula ef(formula);
    return *this - ef;
  }

  EmpiricalFormula & EmpiricalFormula::operator-=(const EmpiricalFormula & formula)
  {
    Map<const Element *, SignedSize>::ConstIterator it = formula.formula_.begin();
    for (; it != formula.formula_.end(); ++it)
    {
      if (formula_.has(it->first))
      {
        formula_[it->first] -= it->second;
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

  EmpiricalFormula & EmpiricalFormula::operator-=(const String & formula)
  {
    Map<const Element *, SignedSize> str_formula;
    SignedSize charge = parseFormula_(str_formula, formula);
    charge_ -= charge;
    Map<const Element *, SignedSize>::ConstIterator it = str_formula.begin();
    for (; it != str_formula.end(); ++it)
    {
      if (formula_.has(it->first))
      {
        formula_[it->first] -= it->second;
      }
      else
      {
        formula_[it->first] = -it->second;
      }
    }
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

  bool EmpiricalFormula::hasElement(const Element * element) const
  {
    return formula_.has(element);
  }

  bool EmpiricalFormula::hasElement(const String & element) const
  {
    if (!element_db_->hasElement(element))
    {
      return false;
    }
    else
    {
      if (formula_.has(element_db_->getElement(element)))
      {
        return true;
      }
    }
    return false;
  }

  bool EmpiricalFormula::hasElement(UInt atomic_number) const
  {
    if (!element_db_->hasElement(atomic_number))
    {
      return false;
    }
    else
    {
      if (formula_.has(element_db_->getElement(atomic_number)))
      {
        return true;
      }
    }
    return false;
  }

  bool EmpiricalFormula::operator==(const EmpiricalFormula & formula) const
  {
    return formula_ == formula.formula_ && charge_ == formula.charge_;
  }

  bool EmpiricalFormula::operator==(const String & formula) const
  {
    Map<const Element *, SignedSize> str_formula;
    SignedSize charge = parseFormula_(str_formula, formula);
    return formula_ == str_formula && charge_ == charge;
  }

  bool EmpiricalFormula::operator!=(const EmpiricalFormula & formula) const
  {
    return formula_ != formula.formula_ || charge_ != formula.charge_;
  }

  bool EmpiricalFormula::operator!=(const String & formula) const
  {
    Map<const Element *, SignedSize> str_formula;
    SignedSize charge = parseFormula_(str_formula, formula);
    return formula_ != str_formula || charge_ != charge;
  }

  ostream & operator<<(ostream & os, const EmpiricalFormula & formula)
  {
    Map<String, SignedSize> new_formula;
    for (Map<const Element *, SignedSize>::ConstIterator it = formula.formula_.begin(); it != formula.formula_.end(); ++it)
    {
      new_formula[it->first->getSymbol()] = it->second;
    }

    for (Map<String, SignedSize>::ConstIterator it = new_formula.begin(); it != new_formula.end(); ++it)
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

  SignedSize EmpiricalFormula::parseFormula_(Map<const Element *, SignedSize> & ef, const String & input_formula) const
  {
    SignedSize charge = 0;
    String formula(input_formula), symbol, number;

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
      String charge_part;
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

      if (element_db_->hasElement(symbol))
      {
        if (num != 0)
        {
          const Element * e = element_db_->getElement(symbol);
          if (ef.has(e))
          {
            ef[e] += num;
          }
          else
          {
            ef[e] = num;
          }
        }
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'" + split + "'", "'" + symbol + "'");
      }
    }

    // remove elements with 0 counts
    Map<const Element *, SignedSize>::iterator it = ef.begin();
    while (it != ef.end())
    {
      if (it->second == 0)
      {
         ef.erase(it++);  // Note: post increment needed! Otherwise iterator is invalidated 
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
    Map<const Element *, SignedSize>::iterator it = formula_.begin();
    while (it != formula_.end())
    {
      if (it->second == 0)
      {
         formula_.erase(it++);  // Note: post increment needed! Otherwise iterator is invalidated 
      }
      else
      {
         ++it;
      }
    }
  }

}
