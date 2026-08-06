// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Nikos Patikos $
// $Authors: Nikos Patikos $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>
#include <fstream>

using namespace std;

namespace OpenMS
{
  IsotopePatternGenerator::IsotopePatternGenerator(double probability_cutoff) :
    min_prob_(probability_cutoff)
  {
  }

  IsotopePatternGenerator::IsotopePatternGenerator() :
    min_prob_(1e-15)
  {
  }
  
  IsotopePatternGenerator::~IsotopePatternGenerator() = default;

  IsotopeDistribution IsotopePatternGenerator::runNeutral(const EmpiricalFormula& neutral_formula) const
  {
    if (neutral_formula.isCharged())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "runNeutral() requires a neutral (uncharged) formula, but '" + neutral_formula.toString() +
        "' has charge " + std::to_string(neutral_formula.getCharge()) + ". For an m/z pattern of an ion, use runIon(neutral_formula, adduct).");
    }
    return run(neutral_formula);
  }

  IsotopeDistribution IsotopePatternGenerator::runIon(const EmpiricalFormula& neutral_formula, const AdductInfo& adduct) const
  {
    // resolve the actual atomic composition of the ion (throws on a charged input
    // or an adduct that removes atoms the molecule does not have)
    const EmpiricalFormula ion = adduct.getIonComposition(neutral_formula);

    // isotope pattern of the (neutral) atomic composition
    IsotopeDistribution result = run(ion);

    // convert neutral mass -> m/z: a positive charge is missing |z| electrons,
    // a negative charge carries |z| extra electrons; then divide by |z|
    const int z = adduct.getCharge();
    const double abs_z = std::abs(z);
    for (auto& peak : result)
    {
      peak.setMZ((peak.getMZ() - z * Constants::ELECTRON_MASS_U) / abs_z);
    }
    return result;
  }

}
