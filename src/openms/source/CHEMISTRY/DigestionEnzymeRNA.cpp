// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser  $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <iostream>

using namespace std;

namespace OpenMS
{
  bool DigestionEnzymeRNA::setValueFromFile(const std::string& key, const std::string& value)
  {
    if (DigestionEnzyme::setValueFromFile(key, value))
    {
      return true;
    }
    if (StringUtils::hasSuffix(key, ":CutsAfter"))
    {
      setCutsAfterRegEx(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":CutsBefore"))
    {
      setCutsBeforeRegEx(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":ThreePrimeGain"))
    {
      setThreePrimeGain(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":FivePrimeGain"))
    {
      setFivePrimeGain(value);
      return true;
    }
    return false;
  }

  void DigestionEnzymeRNA::setCutsAfterRegEx(const std::string& value)
  {
    cuts_after_regex_ = value;
  }

  std::string DigestionEnzymeRNA::getCutsAfterRegEx() const
  {
    return cuts_after_regex_;
  }

  void DigestionEnzymeRNA::setCutsBeforeRegEx(const std::string& value)
  {
    cuts_before_regex_ = value;
  }

  std::string DigestionEnzymeRNA::getCutsBeforeRegEx() const
  {
    return cuts_before_regex_;
  }

  void DigestionEnzymeRNA::setThreePrimeGain(const std::string& value)
  {
    three_prime_gain_ = value;
  }

  std::string DigestionEnzymeRNA::getThreePrimeGain() const
  {
    return three_prime_gain_;
  }

  void DigestionEnzymeRNA::setFivePrimeGain(const std::string& value)
  {
    five_prime_gain_ = value;
  }

  std::string DigestionEnzymeRNA::getFivePrimeGain() const
  {
    return five_prime_gain_;
  }
}
