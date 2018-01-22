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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


#include <boost/math/special_functions/binomial.hpp>

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

  double EmpiricalFormula::calculateTheoreticalIsotopesNumber() const
  {
    double total = 1;
    for(auto const & element : formula_)
    {
      UInt non_trace_isotopes = 0;
      auto const & distr = element.first->getIsotopeDistribution();
      for (auto isotope : distr)
      {
        if (isotope.getIntensity() != 0)
        {
          non_trace_isotopes++;
        }
      }
      if (non_trace_isotopes>1 && element.second!=1)
      {
        total *= boost::math::binomial_coefficient<double>(UInt(element.second), non_trace_isotopes);
      }else
      {
        total *= element.second * non_trace_isotopes;
      }
    }
    return total;
  }

  bool EmpiricalFormula::estimateFromWeightAndCompAndS(double average_weight, UInt S, double C, double H, double N, double O, double P)
  {
    const ElementDB* db = ElementDB::getInstance();

    double remaining_weight = average_weight - S * db->getElement("S")->getAverageWeight();

    // The number of sulfurs is set to 0 because we're explicitly specifying their count.
    // We propagate the return value to let the programmer know if the approximation succeeded
    // without requesting a negative number of hydrogens.
    bool ret = estimateFromWeightAndComp(remaining_weight, C, H, N, O, 0.0, P);

    formula_.at(db->getElement("S")) = S;

    return ret;
  }

  bool EmpiricalFormula::estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P)
  {
    const ElementDB* db = ElementDB::getInstance();

    double avgTotal = (C * db->getElement("C")->getAverageWeight() +
                       H * db->getElement("H")->getAverageWeight() +
                       N * db->getElement("N")->getAverageWeight() +
                       O * db->getElement("O")->getAverageWeight() +
                       S * db->getElement("S")->getAverageWeight() +
                       P * db->getElement("P")->getAverageWeight());

    double factor = average_weight / avgTotal;

    formula_.clear();

    formula_.insert(make_pair(db->getElement("C"), (SignedSize) Math::round(C * factor)));
    formula_.insert(make_pair(db->getElement("N"), (SignedSize) Math::round(N * factor)));
    formula_.insert(make_pair(db->getElement("O"), (SignedSize) Math::round(O * factor)));
    formula_.insert(make_pair(db->getElement("S"), (SignedSize) Math::round(S * factor)));
    formula_.insert(make_pair(db->getElement("P"), (SignedSize) Math::round(P * factor)));

    double remaining_mass = average_weight-getAverageWeight();
    SignedSize adjusted_H = Math::round(remaining_mass / db->getElement("H")->getAverageWeight());

    // It's possible for a very small mass to get a negative value here.
    if (adjusted_H < 0)
    {
      // The approximation can still be useful, but we set the return flag to false to explicitly notify the programmer.
      return false;
    }

    // Only insert hydrogens if their number is not negative.
    formula_.insert(make_pair(db->getElement("H"), adjusted_H));
    // The approximation had no issues.
    return true;
  }

  IsotopeDistribution EmpiricalFormula::getIsotopeDistribution(IsotopePatternGenerator* solver) const
  {
    solver->run(*this);
    return *solver;
  }

  IsotopeDistribution EmpiricalFormula::getConditionalFragmentIsotopeDist(const EmpiricalFormula& precursor, const std::set<UInt>& precursor_isotopes) const
  {
    // A fragment's isotopes can only be as high as the largest isolated precursor isotope.
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end())+1;

    // Treat *this as the fragment molecule
    EmpiricalFormula complementary_fragment = precursor-*this;

    IsotopeDistribution fragment_isotope_dist = getIsotopeDistribution(new CoarseIsotopeDistribution(max_depth));
    IsotopeDistribution comp_fragment_isotope_dist = complementary_fragment.getIsotopeDistribution(new CoarseIsotopeDistribution(max_depth));

    CoarseIsotopeDistribution result;
    result.calcFragmentIsotopeDist(fragment_isotope_dist, comp_fragment_isotope_dist, precursor_isotopes);

    // Renormalize to make these conditional probabilities (conditioned on the isolated precursor isotopes)
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
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, formula, "Cannot parse charge part of formula!");
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
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, formula, "This formula does not begin with an element!");
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
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown element '" + split + "'", "'" + symbol + "' found. Please use only valid element identifiers or modify share/OpenMS/CHEMISTRY/Elements.xml!");
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
