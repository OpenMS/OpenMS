// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Ahmed Khalil $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <boost/math/special_functions/binomial.hpp>

#include <iostream>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  EmpiricalFormula::EmpiricalFormula() :
    charge_(0)
  {}

  EmpiricalFormula::EmpiricalFormula(const String& formula)
  {
    charge_ = parseFormula_(formula_, formula);
  }

  EmpiricalFormula::EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge)
  {
    formula_[element] = number;
    charge_ = charge;
  }

  EmpiricalFormula::~EmpiricalFormula() = default;

  double EmpiricalFormula::getMonoWeight() const
  {
    double weight = Constants::PROTON_MASS_U * charge_;
    for (const auto& it : formula_)
    {
      weight += it.first->getMonoWeight() * (double)it.second;
    }
    return weight;
  }

  double EmpiricalFormula::getLightestIsotopeWeight() const
  {
    double weight = Constants::PROTON_MASS_U * charge_;
    for (const auto& it : formula_)
    {
      // Isotopes should be filled sorted by mz in Elements, so we use
      //  the first element instead of getMin()
      weight += it.first->getIsotopeDistribution()[0].getMZ() * (double)it.second;
    }
    return weight;
  }

  double EmpiricalFormula::getAverageWeight() const
  {
    double weight = Constants::PROTON_MASS_U * charge_;
    for (const auto& it : formula_)
    {
      weight += it.first->getAverageWeight() * (double)it.second;
    }
    return weight;
  }

  double EmpiricalFormula::calculateTheoreticalIsotopesNumber() const
  {
    double total = 1;
    for (const auto& element : formula_)
    {
      UInt non_trace_isotopes = 0;
      const auto& distr = element.first->getIsotopeDistribution();
      for (const auto& isotope : distr)
      {
        if (isotope.getIntensity() != 0)
        {
          non_trace_isotopes++;
        }
      }
      if (non_trace_isotopes>1 && element.second!=1)
      {
        total *= boost::math::binomial_coefficient<double>(UInt(element.second), non_trace_isotopes);
      }
      else
      {
        total *= element.second*non_trace_isotopes;
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

  bool EmpiricalFormula::estimateFromMonoWeightAndComp(double mono_weight, double C, double H, double N, double O, double S, double P)
  {
    const ElementDB* db = ElementDB::getInstance();

    double monoTotal = (C * db->getElement("C")->getMonoWeight() +
                       H * db->getElement("H")->getMonoWeight() +
                       N * db->getElement("N")->getMonoWeight() +
                       O * db->getElement("O")->getMonoWeight() +
                       S * db->getElement("S")->getMonoWeight() +
                       P * db->getElement("P")->getMonoWeight());

    double factor = mono_weight / monoTotal;

    formula_.clear();

    formula_.insert(make_pair(db->getElement("C"), (SignedSize) Math::round(C * factor)));
    formula_.insert(make_pair(db->getElement("N"), (SignedSize) Math::round(N * factor)));
    formula_.insert(make_pair(db->getElement("O"), (SignedSize) Math::round(O * factor)));
    formula_.insert(make_pair(db->getElement("S"), (SignedSize) Math::round(S * factor)));
    formula_.insert(make_pair(db->getElement("P"), (SignedSize) Math::round(P * factor)));

    double remaining_mass = mono_weight-getMonoWeight();
    SignedSize adjusted_H = Math::round(remaining_mass / db->getElement("H")->getMonoWeight());

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

  IsotopeDistribution EmpiricalFormula::getIsotopeDistribution(const IsotopePatternGenerator& solver) const
  {
    return solver.run(*this);
  }

  IsotopeDistribution EmpiricalFormula::getConditionalFragmentIsotopeDist(const EmpiricalFormula& precursor,
                                                                          const std::set<UInt>& precursor_isotopes,
                                                                          const CoarseIsotopePatternGenerator& solver) const
  {
    // A fragment's isotopes can only be as high as the largest isolated precursor isotope.
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end())+1;

    // Treat *this as the fragment molecule
    EmpiricalFormula complementary_fragment = precursor-*this;

    IsotopeDistribution fragment_isotope_dist = getIsotopeDistribution(CoarseIsotopePatternGenerator(max_depth));
    IsotopeDistribution comp_fragment_isotope_dist = complementary_fragment.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_depth));

    IsotopeDistribution result = solver.calcFragmentIsotopeDist(fragment_isotope_dist, comp_fragment_isotope_dist, precursor_isotopes, getMonoWeight());

    // Renormalize to make these conditional probabilities (conditioned on the isolated precursor isotopes)
    result.renormalize();

    return result;
  }

  SignedSize EmpiricalFormula::getNumberOf(const Element* element) const
  {
    const auto& it  = formula_.find(element);
    if (it != formula_.end())
    {
      return it->second;
    }
    return 0;
  }

  SignedSize EmpiricalFormula::getNumberOfAtoms() const
  {
    SignedSize num_atoms(0);
    for (const auto& it : formula_) num_atoms += it.second;
    return num_atoms;
  }

  void EmpiricalFormula::setCharge(Int charge)
  {
    charge_ = charge;
  }

  Int EmpiricalFormula::getCharge() const
  {
    return charge_;
  }

  String EmpiricalFormula::toString() const
  {
    String formula;
    auto formula_map = toMap();
    for (const auto& it : formula_map)
    {
      (formula += it.first) += String(it.second);
    }
    return formula;
  }

  std::map<std::string, int> EmpiricalFormula::toMap() const
  {
    std::map<std::string, int> formula_map;
    for (const auto & it : formula_)
    {
      formula_map[it.first->getSymbol()] = it.second;
    }
    return formula_map;
  }

  EmpiricalFormula EmpiricalFormula::operator*(const SignedSize& times) const
  {
    EmpiricalFormula ef(*this);
    for (const auto& it : formula_) ef.formula_[it.first] *= times;
    ef.charge_ *= times;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula EmpiricalFormula::operator+(const EmpiricalFormula& formula) const
  {
    EmpiricalFormula ef;
    ef.formula_ = formula.formula_;
    for (const auto& it : formula_)
    {
      auto ef_it  = ef.formula_.find(it.first);
      if (ef_it != ef.formula_.end())
      {
        ef_it->second += it.second;
      }
      else
      {
        ef.formula_.insert(it);
      }
    }
    ef.charge_ = charge_ + formula.charge_;
    ef.removeZeroedElements_();
    return ef;
  }

  EmpiricalFormula& EmpiricalFormula::operator+=(const EmpiricalFormula& formula)
  {
    for (const auto& it : formula.formula_)
    {
      auto f_it  = formula_.find(it.first);
      if (f_it != formula_.end())
      {
        f_it->second += it.second;
      }
      else
      {
        formula_.insert(it);
      }
    }
    charge_ += formula.charge_;
    removeZeroedElements_();
    return *this;
  }

  EmpiricalFormula EmpiricalFormula::operator-(const EmpiricalFormula& formula) const
  {
    EmpiricalFormula ef(*this);
    for (const auto& it : formula.formula_)
    {
      const Element* e = it.first;
      SignedSize num = it.second;
      auto ef_it  = ef.formula_.find(e);
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
    for (const auto& it : formula.formula_)
    {
      auto f_it  = formula_.find(it.first);
      if (f_it != formula_.end())
      {
        f_it->second -= it.second;
      }
      else
      {
        formula_[it.first] = -it.second;
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

  bool EmpiricalFormula::contains(const EmpiricalFormula& ef) const
  {
    for (const auto& it : ef)
    {
      if (this->getNumberOf(it.first) < it.second)
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
    for (const auto& it : formula.formula_)
    {
      new_formula[it.first->getSymbol()] = it.second;
    }

    for (const auto& it : new_formula)
    {
      os << it.first;
      if (it.second > 1) os << it.second;
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

  Int EmpiricalFormula::parseFormula_(std::map<const Element*, SignedSize>& ef, const String& input_formula) const
  {
    Int charge{0};
    String formula(input_formula);
    formula.trim();

    // we start with the charge part, read until the begin of the formula or a element symbol occurs
    String suffix;
    for (SignedSize reverse_i(formula.size() - 1); reverse_i >= 0; --reverse_i)
    {
      if (!isalpha(formula[reverse_i]))
      {
        suffix.insert(0,1, formula[reverse_i]); // pre-pend
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

        Int tmp_charge = 1;
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
    std::vector<std::string> splitter;
    splitter.reserve(formula.size() / 2); // reasonable estimate for small formulae like C6H12O6
    if (!formula.empty())
    {
      if (!isdigit(formula[0]) || formula[0] == '(')
      {
        bool is_isotope(false), is_symbol(false);
        bool char_is_upper, is_bracket;
        std::string split;
        for (const auto& curr : formula)
        {
          char_is_upper = isupper(curr);
          is_bracket = (curr == '(');
          if ((char_is_upper && (!is_isotope || is_symbol)) || is_bracket)
          {
            if (!split.empty())
            {
              splitter.push_back(std::move(split));
              is_isotope = false;
              is_symbol = false;
            }
            split = curr;
          }
          else
          {
            split += curr;
          }
          if (is_bracket)
          {
            is_isotope = true;
          }
          if (char_is_upper)
          {
            is_symbol = true;
          }
        }
        splitter.push_back(std::move(split));
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, formula, "This formula does not begin with an element!");
      }
    }

    // add up the elements
    const ElementDB* db = ElementDB::getInstance();
    for (Size i = 0; i != splitter.size(); ++i)
    {
      const String& split = splitter[i];
      String number;
      String symbol;
      bool had_symbol(false);
      for (SignedSize j = split.size() - 1; j >= 0; --j)
      {
        if (!had_symbol && (isdigit(split[j]) || split[j] == '-'))
        {
          number.insert(0,1, split[j]); // pre-pend
        }
        else
        {
          symbol.insert(0,1, split[j]); // pre-pend
          had_symbol = true;
        }
      }

      SignedSize num(1);
      if (!number.empty())
      {
        num = number.toInt();
      }

      const Element* e = db->getElement(symbol);
      if (e != nullptr)
      {
        if (num != 0)
        {
          auto it = ef.find(e);
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
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown element '" + split + "'", "'" + symbol + "' found.");
      }
    }

    // remove elements with 0 counts
    auto it = ef.begin();
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

  bool EmpiricalFormula::operator<(const EmpiricalFormula& rhs) const
  {
    if (formula_.size() != rhs.formula_.size())
    {
      return formula_.size() < rhs.formula_.size();
    }

    if (charge_ != rhs.charge_) 
    {
      return charge_ < rhs.charge_;
    }
    
    return formula_ < rhs.formula_;
  }

  EmpiricalFormula EmpiricalFormula::hydrogen(int n_atoms)
  {
    const ElementDB* db = ElementDB::getInstance();
    return EmpiricalFormula(n_atoms, db->getElement(1));
  }

  EmpiricalFormula EmpiricalFormula::water(int n_molecules)
  {
    const ElementDB* db = ElementDB::getInstance();
    EmpiricalFormula formula;
    formula.formula_[db->getElement(1)] = n_molecules * 2; // hydrogen
    formula.formula_[db->getElement(8)] = n_molecules; // oxygen
    return formula;
  }

} // namespace OpenMS
