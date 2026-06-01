// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  const EmpiricalFormula Ribonucleotide::default_baseloss_ =
    EmpiricalFormula("C5H10O5");

  ostream& operator<<(ostream& os, const Ribonucleotide& ribo)
  {
    os << "Ribonucleotide '"
       << ribo.code_ << "' ("
       << ribo.name_ << ", "
       << ribo.formula_ << ")";
    return os;
  }

  Ribonucleotide::Ribonucleotide(
    const String& name, const String& code, const String& new_code,
    const String& html_code, const EmpiricalFormula& formula, char origin,
    double mono_mass, double avg_mass, enum TermSpecificityNuc term_spec,
    const EmpiricalFormula& baseloss_formula):
    name_(name), code_(code), new_code_(new_code), html_code_(html_code),
    formula_(formula), origin_(origin), mono_mass_(mono_mass),
    avg_mass_(avg_mass), term_spec_(term_spec),
    baseloss_formula_(baseloss_formula)
  {
  }

  Ribonucleotide::~Ribonucleotide() = default;

  bool Ribonucleotide::operator==(const Ribonucleotide& ribonucleotide) const
  {
    return name_ == ribonucleotide.name_ &&
        code_ == ribonucleotide.code_ &&
        new_code_ == ribonucleotide.new_code_ &&
        html_code_ == ribonucleotide.html_code_ &&
        formula_ == ribonucleotide.formula_ &&
        origin_ == ribonucleotide.origin_ &&
        mono_mass_ == ribonucleotide.mono_mass_ &&
        avg_mass_ == ribonucleotide.avg_mass_ &&
        term_spec_ == ribonucleotide.term_spec_ &&
        baseloss_formula_ == ribonucleotide.baseloss_formula_;
  }

  const String Ribonucleotide::getCode() const
  {
    return code_;
  }

  void Ribonucleotide::setCode(const String& code)
  {
    code_ = code;
  }

  const String Ribonucleotide::getName() const
  {
    return name_;
  }

  void Ribonucleotide::setName(const String& name)
  {
    name_ = name;
  }

  double Ribonucleotide::getMonoMass() const
  {
      return mono_mass_;
  }

  void Ribonucleotide::setMonoMass(double mono_mass)
  {
    mono_mass_ = mono_mass;
  }

  double Ribonucleotide::getAvgMass() const
  {
    return avg_mass_;
  }

  void Ribonucleotide::setAvgMass(double avg_mass)
  {
    avg_mass_ = avg_mass;
  }

  const String Ribonucleotide::getNewCode() const
  {
    return new_code_;
  }

  void Ribonucleotide::setNewCode(const String& new_code)
  {
    new_code_ = new_code;
  }

  char Ribonucleotide::getOrigin() const
  {
    return origin_;
  }

  void Ribonucleotide::setOrigin(char origin)
  {
    origin_ = origin;
  }

  String Ribonucleotide::getHTMLCode() const
  {
    return html_code_;
  }

  void Ribonucleotide::setHTMLCode(const String& html_code)
  {
    html_code_ = html_code;
  }

  const EmpiricalFormula Ribonucleotide::getFormula() const
  {
    return formula_;
  }

  void Ribonucleotide::setFormula(const EmpiricalFormula& formula)
  {
    formula_ = formula;
  }

  enum Ribonucleotide::TermSpecificityNuc Ribonucleotide::getTermSpecificity() const
  {
    return term_spec_;
  }

  void Ribonucleotide::setTermSpecificity(enum TermSpecificityNuc term_spec)
  {
    if (term_spec == NUMBER_OF_TERM_SPECIFICITY)
    {
      String msg = "invalid terminal specificity";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg, "NUMBER_OF_TERM_SPECIFICITY");
    }
    term_spec_ = term_spec;
  }

  const EmpiricalFormula Ribonucleotide::getBaselossFormula() const
  {
    return baseloss_formula_;
  }

  void Ribonucleotide::setBaselossFormula(const EmpiricalFormula& formula)
  {
    baseloss_formula_ = formula;
  }

  bool Ribonucleotide::isModified() const
  {
    return (code_.length() != 1) || (code_[0] != origin_);
  }

  bool Ribonucleotide::isAmbiguous() const
  {
    return code_.back() == '?';
  }

}
