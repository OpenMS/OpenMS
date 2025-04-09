// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Adduct.h>

#include <OpenMS/CHEMISTRY/Element.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <iostream>

#include <map>

namespace OpenMS
{

  Adduct::Adduct() :
    charge_(0),
    amount_(0),
    singleMass_(0),
    log_prob_(0),
    formula_(),
    rt_shift_(0),
    label_()
  {
  }

  Adduct::Adduct(Int charge) :
    charge_(charge),
    amount_(0),
    singleMass_(0),
    log_prob_(0),
    formula_(),
    rt_shift_(0),
    label_()
  {
  }

  Adduct::Adduct(Int charge, Int amount, double singleMass, const String& formula, double log_prob, double rt_shift, const String& label) :
    charge_(charge),
    amount_(amount),
    singleMass_(singleMass),
    log_prob_(log_prob),
    rt_shift_(rt_shift),
    label_(label)
  {
    if (amount < 0)
    {
      std::cerr << "Attention: Adduct received negative amount! (" << amount << ")\n";
    }
    formula_ = checkFormula_(formula);
  }

  Adduct Adduct::operator*(const Int m) const
  {
    Adduct a = *this;
    a.amount_ *= m;
    return a;
  }

  Adduct Adduct::operator+(const Adduct& rhs)
  {
    if (this->formula_ != rhs.formula_)
    {
      throw "Adduct::Operator +()  tried to add incompatible adduct!";
    }
    Adduct a = *this;
    a.amount_ += rhs.amount_;
    return a;
  }

  void Adduct::operator+=(const Adduct& rhs)
  {
    if (this->formula_ != rhs.formula_)
    {
      throw "Adduct::Operator +=()  tried to add incompatible adduct!";
    }
    this->amount_ += rhs.amount_;
  }

  //@{ Accessors
  const Int& Adduct::getCharge() const
  {
    return charge_;
  }

  void Adduct::setCharge(const Int& charge)
  {
    charge_ = charge;
  }

  const Int& Adduct::getAmount() const
  {
    return amount_;
  }

  void Adduct::setAmount(const Int& amount)
  {
    if (amount < 0)
    {
      std::cerr << "Warning: Adduct received negative amount! (" << amount << ")\n";
    }
    amount_ = amount;
  }

  const double& Adduct::getSingleMass() const
  {
    return singleMass_;
  }

  void Adduct::setSingleMass(const double& singleMass)
  {
    singleMass_ = singleMass;
  }

  const double& Adduct::getLogProb() const
  {
    return log_prob_;
  }

  void Adduct::setLogProb(const double& log_prob)
  {
    log_prob_ = log_prob;
  }

  const String& Adduct::getFormula() const
  {
    return formula_;
  }

  void Adduct::setFormula(const String& formula)
  {
    formula_ = checkFormula_(formula);
  }

  const double& Adduct::getRTShift() const
  {
    return rt_shift_;
  }

  const String& Adduct::getLabel() const
  {
    return label_;
  }
  
  String Adduct::toAdductString(const String& ion_string, const Int& charge)
  {
    EmpiricalFormula ef(ion_string);
    String charge_sign = charge >= 0 ? "+" : "-";
    String s("[M");

    //need elements sorted canonically (by string)
    std::map<String, String> sorted_elem_map;
    for (const auto& element_count : ef)
    {
      String e_symbol(element_count.first->getSymbol());
      String tmp = element_count.second > 0 ? "+" : "-";
      tmp += std::abs(element_count.second) > 1 ? String(std::abs(element_count.second)) : "";
      tmp += e_symbol;
      sorted_elem_map[e_symbol] = std::move(tmp);
    }
    for (const auto& sorted_e_cnt : sorted_elem_map)
    {
      s += sorted_e_cnt.second;
    }
    s += String("]");
    s += std::abs(charge) > 1 ? String(std::abs(charge)) : "";
    s += charge_sign;

    return s;
  }

  String Adduct::checkFormula_(const String& formula)
  {
    EmpiricalFormula ef(formula);
    if (ef.getCharge() != 0)
    {
      std::cerr << "Warning: Adduct contains explicit charge (alternating mass)! (" << formula << ")\n";
    }
    if (ef.isEmpty())
    {
      std::cerr << "Warning: Adduct was given empty formula! (" << formula << ")\n";
    }
    if ((ef.getNumberOfAtoms() > 1) && (std::distance(ef.begin(), ef.end()) == 1))
    {
      std::cerr << "Warning: Adduct was given only a single element but with an abundance>1. This might lead to errors! (" << formula << ")\n";
    }

    return ef.toString();
  }

  ///Print the contents of an Adduct to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Adduct& a)
  {
    os << "---------- Adduct -----------------\n";
    os << "Charge: " << a.charge_ << std::endl;
    os << "Amount: " << a.amount_ << std::endl;
    os << "MassSingle: " << a.singleMass_ << std::endl;
    os << "Formula: " << a.formula_ << std::endl;
    os << "log P: " << a.log_prob_ << std::endl;
    return os;
  }

  OPENMS_DLLAPI bool operator==(const Adduct& a, const  Adduct& b)
  {
    return a.charge_ == b.charge_
           && a.amount_ == b.amount_
           && a.singleMass_ == b.singleMass_
           && a.log_prob_ == b.log_prob_
           && a.formula_ == b.formula_;

  }

}
