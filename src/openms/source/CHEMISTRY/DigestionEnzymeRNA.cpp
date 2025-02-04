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
  bool DigestionEnzymeRNA::setValueFromFile(const String& key, const String& value)
  {
    if (DigestionEnzyme::setValueFromFile(key, value))
    {
      return true;
    }
    if (key.hasSuffix(":CutsAfter"))
    {
      setCutsAfterRegEx(value);
      return true;
    }
    if (key.hasSuffix(":CutsBefore"))
    {
      setCutsBeforeRegEx(value);
      return true;
    }
    if (key.hasSuffix(":ThreePrimeGain"))
    {
      setThreePrimeGain(value);
      return true;
    }
    if (key.hasSuffix(":FivePrimeGain"))
    {
      setFivePrimeGain(value);
      return true;
    }
    return false;
  }

  void DigestionEnzymeRNA::setCutsAfterRegEx(const String& value)
  {
    cuts_after_regex_ = value;
  }

  String DigestionEnzymeRNA::getCutsAfterRegEx() const
  {
    return cuts_after_regex_;
  }

  void DigestionEnzymeRNA::setCutsBeforeRegEx(const String& value)
  {
    cuts_before_regex_ = value;
  }

  String DigestionEnzymeRNA::getCutsBeforeRegEx() const
  {
    return cuts_before_regex_;
  }

  void DigestionEnzymeRNA::setThreePrimeGain(const String& value)
  {
    three_prime_gain_ = value;
  }

  String DigestionEnzymeRNA::getThreePrimeGain() const
  {
    return three_prime_gain_;
  }

  void DigestionEnzymeRNA::setFivePrimeGain(const String& value)
  {
    five_prime_gain_ = value;
  }

  String DigestionEnzymeRNA::getFivePrimeGain() const
  {
    return five_prime_gain_;
  }
}
