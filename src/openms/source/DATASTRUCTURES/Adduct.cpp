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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <iostream>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

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

  Adduct::Adduct(Int charge, Int amount, double singleMass, String formula, double log_prob, double rt_shift, const String label) :
    charge_(charge),
    amount_(amount),
    singleMass_(singleMass),
    log_prob_(log_prob),
    rt_shift_(rt_shift),
    label_(label)
  {
    if (amount < 0)
      std::cerr << "Attention: Adduct received negative amount! (" << amount << ")\n";
    formula_ = checkFormula_(formula);
  }

  Adduct Adduct::operator*(const Int m) const
  {
    Adduct a = *this;
    a.amount_ *= m;
    return a;
  }

  Adduct Adduct::operator+(const Adduct & rhs)
  {
    if (this->formula_ != rhs.formula_)
    {
      throw "Adduct::Operator +()  tried to add incompatible adduct!";
    }
    Adduct a = *this;
    a.amount_ += rhs.amount_;
    return a;
  }

  void Adduct::operator+=(const Adduct & rhs)
  {
    if (this->formula_ != rhs.formula_)
    {
      throw "Adduct::Operator +=()  tried to add incompatible adduct!";
    }
    this->amount_ += rhs.amount_;
  }

  //@{ Accessors
  const Int & Adduct::getCharge() const
  {
    return charge_;
  }

  void Adduct::setCharge(const Int & charge)
  {
    charge_ = charge;
  }

  const Int & Adduct::getAmount() const
  {
    return amount_;
  }

  void Adduct::setAmount(const Int & amount)
  {
    if (amount < 0)
      std::cerr << "Warning: Adduct received negative amount! (" << amount << ")\n";
    amount_ = amount;
  }

  const double & Adduct::getSingleMass() const
  {
    return singleMass_;
  }

  void Adduct::setSingleMass(const double & singleMass)
  {
    singleMass_ = singleMass;
  }

  const double & Adduct::getLogProb() const
  {
    return log_prob_;
  }

  void Adduct::setLogProb(const double & log_prob)
  {
    log_prob_ = log_prob;
  }

  const String & Adduct::getFormula() const
  {
    return formula_;
  }

  void Adduct::setFormula(const String & formula)
  {
    formula_ = checkFormula_(formula);
  }

  const double & Adduct::getRTShift() const
  {
    return rt_shift_;
  }

  const String & Adduct::getLabel() const
  {
    return label_;
  }

  String Adduct::checkFormula_(const String & formula)
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
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Adduct & a)
  {
    os  << "---------- Adduct -----------------\n";
    os << "Charge: " << a.charge_ << std::endl;
    os << "Amount: " << a.amount_ << std::endl;
    os << "MassSingle: " << a.singleMass_ << std::endl;
    os << "Formula: " << a.formula_ << std::endl;
    os << "log P: " << a.log_prob_ << std::endl;
    return os;
  }

  OPENMS_DLLAPI bool operator==(const Adduct & a, const  Adduct & b)
  {
    return a.charge_ == b.charge_
           && a.amount_ == b.amount_
           && a.singleMass_ == b.singleMass_
           && a.log_prob_ == b.log_prob_
           && a.formula_ == b.formula_;

  }

}
