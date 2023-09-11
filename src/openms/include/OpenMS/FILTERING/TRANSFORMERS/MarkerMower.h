// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
    @brief MarkerMower uses PeakMarker to find peaks, those that are not marked get removed

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI MarkerMower :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    MarkerMower();
    /// destructor
    ~MarkerMower() override;

    /// copy constructor
    MarkerMower(const MarkerMower & source);
    /// assignment operator
    MarkerMower & operator=(const MarkerMower & source);
    // @}

    // @name Accessors
    // @{
    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      typedef typename SpectrumType::Iterator Iterator;

      std::map<double, int> marks;
      for (std::vector<PeakMarker *>::const_iterator cvit = markers_.begin(); cvit != markers_.end(); ++cvit)
      {
        std::map<double, bool> marked;
        (*cvit)->apply(marked, spectrum);
        for (std::map<double, bool>::const_iterator cmit = marked.begin(); cmit != marked.end(); ++cmit)
        {
          if (cmit->second)
          {
            marks[cmit->first]++;
          }
        }
      }

      for (Iterator it = spectrum.begin(); it != spectrum.end(); )
      {
        if (marks[it->getMZ()] > 0)
        {
          ++it;
        }
        else
        {
          it = spectrum.erase(it);
        }
      }
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    static const String getProductName()
    {
      return "MarkerMower";
    }

    /// insert new Marker (violates the DefaultParamHandler interface)
    void insertmarker(PeakMarker * peak_marker);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    // @}

private:
    /// used peak markers
    std::vector<PeakMarker *> markers_;

  };
}
