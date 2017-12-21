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

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerSH.h>

#include <vector>

namespace OpenMS
{
  PeakPickerSH::PeakPickerSH() :
    DefaultParamHandler("PeakPickerSH"),
    ProgressLogger()
  {
    defaultsToParam_();
  }

  PeakPickerSH::~PeakPickerSH()
  {
    // FLO: Do not care at the moment
  }

  void PeakPickerSH::pickExperiment(const PeakMap & input, PeakMap & output)
  {
    // make sure that output is clear
    output.clear(true);

    // copy experimental settings
    static_cast<ExperimentalSettings &>(output) = input;

    // resize output with respect to input
    output.resize(input.size());

    std::cout << "Before loop, input size = " << input.size() << std::endl;
    Size progress = 0;
    for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
    {
      output[scan_idx].clear(true);
      output[scan_idx].SpectrumSettings::operator=(input[scan_idx]);
      output[scan_idx].MetaInfoInterface::operator=(input[scan_idx]);
      output[scan_idx].setRT(input[scan_idx].getRT());
      output[scan_idx].setMSLevel(input[scan_idx].getMSLevel());
      output[scan_idx].setName(input[scan_idx].getName());
      output[scan_idx].setType(SpectrumSettings::CENTROID);

      if (input[scan_idx].getMSLevel() != 1)
      {
        // When not considering MS2 data (MS2 fragment mass tracing=0), Lukas leaves out
        // the entire scan (instead of just copying it to the output as seen in
        // another plugin).
        // pick(input[scan_idx], output[scan_idx], 4.0);
      }
      else
      {
        // TODO: Read value 4.0 from parameters # PeakPickerSH.cpp
        pick(input[scan_idx], output[scan_idx], 5.0);
      }
      setProgress(++progress);
    }
    std::cout << "After loop" << std::endl;

    endProgress();
  }

}
