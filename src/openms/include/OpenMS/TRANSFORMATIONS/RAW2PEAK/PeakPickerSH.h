// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERSH_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERSH_H

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

#endif
