// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_DATAACCESSHELPER_H_
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_DATAACCESSHELPER_H_

#include <fstream>

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenMS
{
  /**
    @brief Several helpers to convert OpenMS datastructures to structures that
           implement the OpenSWATH interfaces.
  */
  class OPENMS_DLLAPI OpenSwathDataAccessHelper
  {
public:
    /// Convert a SpectrumPtr to an OpenMS Spectrum
    static void convertToOpenMSSpectrum(OpenMS::MSSpectrum<> & spectrum, const OpenSwath::SpectrumPtr sptr);

    static OpenSwath::SpectrumPtr convertToSpectrumPtr(const OpenMS::MSSpectrum<> & spectrum)
    {
      // const MSSpectrumType & spectrum = (*ms_experiment_)[id];
      OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
      OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
      for (MSSpectrum<>::const_iterator it = spectrum.begin(); it != spectrum.end(); it++)
      {
        mz_array->data.push_back(it->getMZ());
        intensity_array->data.push_back(it->getIntensity());
      }

      // push back mz first, then intensity.
      // FEATURE (hroest) annotate which is which
      //std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
      //binaryDataArrayPtrs.push_back(mz_array);
      //binaryDataArrayPtrs.push_back(intensity_array);

      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      //sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;
      sptr->setMZArray(mz_array);
      sptr->setIntensityArray(intensity_array);

      return sptr;
    }

    /// Convert a ChromatogramPtr to an OpenMS Chromatogram
    static void convertToOpenMSChromatogram(OpenMS::MSChromatogram<> & chromatogram, const OpenSwath::ChromatogramPtr cptr);

    /// convert from the OpenMS Targeted experiment to the light Targeted Experiment
    static void convertTargetedExp(OpenMS::TargetedExperiment & transition_exp_, OpenSwath::LightTargetedExperiment & transition_exp);

  };

} //end namespace OpenMS

#endif
