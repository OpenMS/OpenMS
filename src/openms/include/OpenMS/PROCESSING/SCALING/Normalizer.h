// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Normalizes the peak intensities spectrum-wise.

    Either to a total intensity-sum of one (i.e. to total-ion-count; TIC) or to a maximum intensity of one.

    @htmlinclude OpenMS_Normalizer.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI Normalizer :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    Normalizer();
    /// destructor
    ~Normalizer() override;

    /// assignment operator
    Normalizer & operator=(const Normalizer & source);
    /// copy constructor
    Normalizer(const Normalizer & source);

    // @}

    // @name Accessors
    // @{

    /**
      @brief Workhorse of this class.

      @param spectrum Input/output spectrum containing peaks
      @throws Exception::InvalidValue if 'method_' has unknown value
    */
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType& spectrum) const
    {
      if (spectrum.empty()) return;

      typedef typename SpectrumType::Iterator Iterator;
      typedef typename SpectrumType::ConstIterator ConstIterator;

      double divisor(0);
      // find divisor      
      if (method_ == "to_one")
      { // normalizes the max peak to 1 and the remaining peaks to values relative to max
        divisor = spectrum.begin()->getIntensity(); // safety measure: if all intensities are negative, divisor would stay 0 (as constructed)
        for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
        {
          if (divisor < it->getIntensity()) divisor = it->getIntensity();
        }
      }
      else if (method_ == "to_TIC")
      { // normalizes the peak intensities to the TIC
        for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
        {
          divisor += it->getIntensity();
        }
      }
      // method unknown
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Method not known", method_);
      }

      // normalize
      for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        it->setIntensity(it->getIntensity() / divisor);
      }

      return;

    }

    ///
    void filterPeakSpectrum(PeakSpectrum & spectrum) const;
    ///
    void filterPeakMap(PeakMap & exp) const;

    void updateMembers_() override;

    // @}

private:
    String method_;
  };


}
