// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  DigestionEnzymeProtein::DigestionEnzymeProtein() :
    DigestionEnzyme(),
    n_term_gain_(""),
    c_term_gain_(""),
    psi_id_(""),
    xtandem_id_(""),
    comet_id_(-1),
    msgf_id_(-1),
    omssa_id_(-1)
  {
  }

  DigestionEnzymeProtein::DigestionEnzymeProtein(const DigestionEnzyme& d) :
      DigestionEnzyme(d),
      n_term_gain_(""),
      c_term_gain_(""),
      psi_id_(""),
      xtandem_id_(""),
      comet_id_(-1),
      msgf_id_(-1),
      omssa_id_(-1)
  {
  }

  DigestionEnzymeProtein::DigestionEnzymeProtein(const std::string& name,
                                                 const std::string& cleavage_regex,
                                                 const std::set<std::string>& synonyms,
                                                 std::string regex_description,
                                                 EmpiricalFormula n_term_gain,
                                                 EmpiricalFormula c_term_gain,
                                                 std::string psi_id,
                                                 std::string xtandem_id,
                                                 Int comet_id,
                                                 Int msgf_id,
                                                 Int omssa_id) :
    DigestionEnzyme(name, cleavage_regex, synonyms, std::move(regex_description)),
    n_term_gain_(std::move(n_term_gain)),
    c_term_gain_(std::move(c_term_gain)),
    psi_id_(std::move(psi_id)),
    xtandem_id_(std::move(xtandem_id)),
    comet_id_(comet_id),
    msgf_id_(msgf_id),
    omssa_id_(omssa_id)
  {
  }

  DigestionEnzymeProtein::~DigestionEnzymeProtein() = default;

  void DigestionEnzymeProtein::setNTermGain(const EmpiricalFormula& value)
  {
    n_term_gain_ = value;
  }

  void DigestionEnzymeProtein::setCTermGain(const EmpiricalFormula& value)
  {
    c_term_gain_ = value;
  }

  EmpiricalFormula DigestionEnzymeProtein::getNTermGain() const
  {
    return n_term_gain_;
  }

  EmpiricalFormula DigestionEnzymeProtein::getCTermGain() const
  {
    return c_term_gain_;
  }

  void DigestionEnzymeProtein::setPSIID(const std::string& value)
  {
    psi_id_ = value;
  }

  std::string DigestionEnzymeProtein::getPSIID() const
  {
    return psi_id_;
  }

  void DigestionEnzymeProtein::setXTandemID(const std::string& value)
  {
    xtandem_id_ = value;
  }

  std::string DigestionEnzymeProtein::getXTandemID() const
  {
    return xtandem_id_;
  }

  void DigestionEnzymeProtein::setCometID(Int value)
  {
    comet_id_ = value;
  }

  Int DigestionEnzymeProtein::getCometID() const
  {
    return comet_id_;
  }

  void DigestionEnzymeProtein::setOMSSAID(Int value)
  {
    omssa_id_ = value;
  }

  Int DigestionEnzymeProtein::getOMSSAID() const
  {
    return omssa_id_;
  }

  void DigestionEnzymeProtein::setMSGFID(Int value)
  {
    msgf_id_ = value;
  }

  Int DigestionEnzymeProtein::getMSGFID() const
  {
    return msgf_id_;
  }

  bool DigestionEnzymeProtein::operator==(const DigestionEnzymeProtein& enzyme) const
  {
    return DigestionEnzyme::operator==(enzyme) &&
           n_term_gain_ == enzyme.n_term_gain_ &&
           c_term_gain_ == enzyme.c_term_gain_ &&
           psi_id_ == enzyme.psi_id_ &&
           xtandem_id_ == enzyme.xtandem_id_ &&
           comet_id_ == enzyme.comet_id_ &&
           msgf_id_ == enzyme.msgf_id_ &&
           omssa_id_ == enzyme.omssa_id_;
  }

  // Note: comparison operators are not inherited. TODO rename it and make virtual
  bool DigestionEnzymeProtein::operator==(const std::string& cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool DigestionEnzymeProtein::operator!=(const std::string& cleavage_regex) const
  {
    return cleavage_regex_ != cleavage_regex;
  }

  bool DigestionEnzymeProtein::operator!=(const DigestionEnzymeProtein& enzyme) const
  {
    return !(*this == enzyme);
  }

  bool DigestionEnzymeProtein::operator<(const DigestionEnzymeProtein& enzyme) const
  {
    return this->getName() < enzyme.getName();
  }

  bool DigestionEnzymeProtein::setValueFromFile(const std::string& key, const std::string& value)
  {
    if (DigestionEnzyme::setValueFromFile(key, value))
    {
      return true;
    }
    if (StringUtils::hasSuffix(key, ":NTermGain"))
    {
      setNTermGain(EmpiricalFormula(value));
      return true;
    }
    if (StringUtils::hasSuffix(key, ":CTermGain"))
    {
      setCTermGain(EmpiricalFormula(value));
      return true;
    }
    if (StringUtils::hasSuffix(key, ":PSIID"))
    {
      setPSIID(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":XTandemID"))
    {
      setXTandemID(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":CometID"))
    {
      setCometID(StringUtils::toInt32(value));
      return true;
    }
    if (StringUtils::hasSuffix(key, ":OMSSAID"))
    {
      setOMSSAID(StringUtils::toInt32(value));
      return true;
    }
    if (StringUtils::hasSuffix(key, ":MSGFID"))
    {
      setMSGFID(StringUtils::toInt32(value));
      return true;
    }
    return false;
  }

  ostream& operator<<(ostream& os, const DigestionEnzymeProtein& enzyme)
  {
    os << static_cast<const DigestionEnzyme&>(enzyme) << " "
       << enzyme.psi_id_;
    return os;
  }

}

