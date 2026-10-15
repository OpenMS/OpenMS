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
#include <OpenMS/CONCEPT/LogStream.h>

#include <atomic>

namespace OpenMS
{
  IsotopeDistribution FineIsotopePatternGenerator::run(const EmpiricalFormula& formula) const
  {
    if (formula.getCharge() < 0)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "FineIsotopePatternGenerator does not support negative charges (formula: " + formula.toString() + ").");
    }
    if (formula.getCharge() > 0)
    {
      // DEPRECATED (OpenMS 3.x): a non-zero charge implicitly adds 'charge'-many H atoms to shift the pattern.
      // This will change in OpenMS 4.0, where the charge is ignored and the neutral pattern is returned.
      // run() is const and may be called from parallel regions, so the warn-once flag must be atomic.
      static std::atomic<bool> warned_once{false};
      if (!warned_once.exchange(true))
      {
        OPENMS_LOG_WARN << "Warning: FineIsotopePatternGenerator::run() was called with a non-zero charge ("
                        << formula.getCharge() << "). The generator currently adds 'charge'-many hydrogen atoms to "
                        << "shift the isotope pattern. This is deprecated and will change in OpenMS 4.0, where the "
                        << "charge will be ignored and the neutral pattern returned. To keep the current behavior, make "
                        << "the adduct explicit via EmpiricalFormula::addChargeAdduct(charge) and run() on the resulting "
                        << "(neutral) formula. To obtain the neutral pattern now, set the charge to 0. "
                        << "(This warning is shown once.)" << std::endl;
      }

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

