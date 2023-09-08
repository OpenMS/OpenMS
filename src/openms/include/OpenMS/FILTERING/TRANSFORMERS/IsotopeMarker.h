// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <map>
#include <vector>
#include <cmath>
#include <utility>

namespace OpenMS
{

  /**
      @brief IsotopeMarker marks peak pairs which could represent an ion and its isotope

        @todo implement a real isotope marking here with isotopedistributions and fitting (Andreas)

        @htmlinclude OpenMS_IsotopeMarker.parameters

        @ingroup PeakMarker
  */
  class OPENMS_DLLAPI IsotopeMarker :
    public PeakMarker
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    IsotopeMarker();

    /// copy constructor
    IsotopeMarker(const IsotopeMarker & source);

    /// destructor
    ~IsotopeMarker() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    IsotopeMarker & operator=(const IsotopeMarker & source);
    // @}

    // @name Accessors
    // @{
    ///
    static PeakMarker * create() { return new IsotopeMarker(); }

    ///
    template <typename SpectrumType>
    void apply(std::map<double, bool> & marked, SpectrumType & spectrum)
    {
      double mzvariation = (double)param_.getValue("mz_variation");
      double invariation = (double)param_.getValue("in_variation");
      Size marks = param_.getValue("marks");

      spectrum.sortByPosition();

      std::map<double, Size> isotopemarks;        // possible isotopes

      for (Size i = 0; i < spectrum.size(); ++i)
      {
        double mz = spectrum[i].getPosition()[0];
        double intensity = spectrum[i].getIntensity();
        Size j = i + 1;

        //std::vector<std::pair<double, double> > isotopes = SpectrumGenerator::instance()->isotopepeaks(mz, intensity);
        CoarseIsotopePatternGenerator solver;
        auto id = solver.estimateFromPeptideWeight(mz);

        while (j < spectrum.size() && spectrum[j].getPosition()[0] <= mz + 3 + mzvariation)
        {
          double curmz = spectrum[j].getPosition()[0];
          double curIntensity = spectrum[j].getIntensity();
          UInt iso = (UInt)(curmz - mz + 0.499999);
          if (iso > 0 && curmz - mz - iso > mzvariation)
          {
            ++j;
            continue;
          }
          if (std::fabs(id.begin()->getIntensity() * intensity - curIntensity) < invariation * id.begin()->getIntensity() * intensity)
          {
            isotopemarks[mz]++;
            isotopemarks[curmz]++;
          }
          ++j;
        }
      }

      for (std::map<double, Size>::const_iterator cmit = isotopemarks.begin(); cmit != isotopemarks.end(); ++cmit)
      {
        if (cmit->second >= marks)
        {
          marked.insert(std::pair<double, bool>(cmit->first, true));
        }
      }
      return;
    }

    ///
    static const String getProductName()
    {
      return "IsotopeMarker";
    }

    // @}

  };

}

