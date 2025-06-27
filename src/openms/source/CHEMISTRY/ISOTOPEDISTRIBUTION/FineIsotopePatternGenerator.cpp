// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>

namespace OpenMS
{
  IsotopeDistribution FineIsotopePatternGenerator::run(const EmpiricalFormula& formula) const
  {
    if (formula.getCharge() != 0)
    {
      // add hydrogen atoms to the formula to match the charge
      EmpiricalFormula charged_formula = formula;
      charged_formula += EmpiricalFormula(formula.getCharge(), ElementDB::getInstance()->getElement("H"));
      charged_formula.setCharge(0); // reset charge, since we added H atoms to match the charge
      /// note: technically, the masses are off by q*electron mass (do we care?)
      return run(charged_formula);
    }

    if (use_total_prob_)
    {
        IsotopeDistribution result(IsoSpecTotalProbWrapper(formula, 1.0-stop_condition_, true).run());
        result.sortByMass();
        return result;
    }
    else
    {
        IsotopeDistribution result(IsoSpecThresholdWrapper(formula, stop_condition_, absolute_).run());
        result.sortByMass();
        return result;
    }
  }

}

