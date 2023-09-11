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

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief ComplementMarker marks peak pairs which could represent y - b ion pairs

        @htmlinclude OpenMS_ComplementMarker.parameters

        @ingroup PeakMarker
  */
  class OPENMS_DLLAPI ComplementMarker :
    public PeakMarker
  {
public:

    // @name Constructors and Destructors
    //@{
    /// standard constructor
    ComplementMarker();

    /// copy constructor
    ComplementMarker(const ComplementMarker & source);

    /// destructor
    ~ComplementMarker() override;
    //@}

    // @name Operators
    //@{
    /// assignment operator
    ComplementMarker & operator=(const ComplementMarker & source);
    //@}

    // @name Accessors
    //@{
    ///
    static PeakMarker * create() { return new ComplementMarker(); }

    ///
    template <typename SpectrumType>
    void apply(std::map<double, bool> marked, SpectrumType & spectrum)
    {
      if (spectrum.size() < 2)
      {
        return;
      }

      // how often a peak needs to be marked to be returned
      double marks = (double)param_.getValue("marks");
      double parentmass = 0.0;
      if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
      double tolerance = (double)param_.getValue("tolerance");
      std::map<double, int> matching_b_y_ions;

      spectrum.sortByPosition();

      SignedSize j = spectrum.size() - 1;
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        while (j >= 0 && spectrum[j].getPosition()[0] > (parentmass - spectrum[i].getPosition()[0]) + tolerance)
        {
          j--;
        }

        // just takes the first matching ion; todo take all
        if (j >= 0 && std::fabs(spectrum[i].getPosition()[0] + spectrum[j].getPosition()[0] - parentmass) < tolerance)
        {
          matching_b_y_ions[spectrum[i].getPosition()[0]]++;
          matching_b_y_ions[spectrum[j].getPosition()[0]]++;
          j--;
        }
      }

      for (std::map<double, int>::const_iterator cmit = matching_b_y_ions.begin(); cmit != matching_b_y_ions.end(); ++cmit)
      {
        if (cmit->second >= marks)
        {
          marked.insert(std::pair<double, bool>(cmit->first, true));
        }
      }
    }

    /// returns the name to register at the factory
    static const String getProductName()
    {
      return "ComplementMarker";
    }

    //@}

  };

}
