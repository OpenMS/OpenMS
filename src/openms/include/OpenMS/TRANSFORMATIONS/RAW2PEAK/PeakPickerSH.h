// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

// TODO: Check if I need this # PeakPickerSH.h
#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING

namespace OpenMS
{
  class OPENMS_DLLAPI PeakPickerSH :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    PeakPickerSH();

    ~PeakPickerSH() override;

    /**
     @brief Picks peaks in one spectrum.
    */
    void pick(const MSSpectrum & input, MSSpectrum & output, float fWindowWidth)
    {
      int i, hw, j;
      double cm, toti, min_dh;

      // Hack: Prepare data structures for Lukas' algorithm
      std::vector<double> masses, intens;
      // TODO: Probably we could save some time when we resize the vectors... # PeakPickerSH.cpp
      //masses.resize(input.size());
      //intens.resize(input.size());
      for (Size k = 0; k < input.size() - 1; ++k)
      {
        // Lukas requires a minimum of intensity (=50). His vectors do not contain
        // other data, so I strip the low ones out right here.
        // TODO: Read 50.0 from parameters  # PeakPickerSH.cpp
        if (input[k].getIntensity() >= 50.0)
        {
          masses.push_back(input[k].getMZ());
          intens.push_back(input[k].getIntensity());
        }
      }

      min_dh = 50.0;              // min height
      hw = fWindowWidth / 2;

      for (i = 2; i < (int)masses.size() - 2; i++)
      {

        // Peak must be concave in the interval [i-2 .. i+2]
        if (intens[i] > min_dh && intens[i] > intens[i - 1] + min_dh && intens[i] >= intens[i + 1] && intens[i - 1] > intens[i - 2] + min_dh && intens[i + 1] >= intens[i + 2])
        {

          cm = 0.0;                   // centroid mass:
          toti = 0.0;             // total intensity:

          for (j = -hw; j <= hw; j++)
          {
            double inte = intens[i - j];
            double mz = masses[i - j];

            cm += inte * mz;
            toti += (double) intens[i - j];
          }
          cm = cm / toti;           // Centre of gravity = centroid

          Peak1D peak;
          peak.setMZ(cm);
          peak.setIntensity(intens[i]);
          output.push_back(peak);
        }
      }
    }

    /**
     @brief Applies the peak-picking algorithm to a map (MSExperiment).

     This method picks peaks for each scan in the map consecutively. The
     resulting picked peaks are written to the output map.
    */
    void pickExperiment(const PeakMap & input, PeakMap & output);
  };
}

