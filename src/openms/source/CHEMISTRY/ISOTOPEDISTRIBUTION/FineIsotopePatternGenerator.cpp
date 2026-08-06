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
      // This will change in OpenMS 4.0, where run() will reject charged formulas (use runNeutral()/runIon()).
      // run() is const and may be called from parallel regions, so the warn-once flag must be atomic.
      static std::atomic<bool> warned_once{false};
      if (!warned_once.exchange(true))
      {
        OPENMS_LOG_WARN << "Warning: FineIsotopePatternGenerator::run() was called with a positively charged "
                        << "EmpiricalFormula (charge " << formula.getCharge() << "). In OpenMS 3.x this implicitly "
                        << "adds 'charge'-many hydrogen atoms and returns a shifted mass pattern (not m/z). This is "
                        << "deprecated; in OpenMS 4.0 run() will reject charged formulas. Use "
                        << "runNeutral(neutral_formula) for a neutral-mass pattern, or "
                        << "runIon(neutral_formula, adduct) for an m/z pattern with explicit ionization. "
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

