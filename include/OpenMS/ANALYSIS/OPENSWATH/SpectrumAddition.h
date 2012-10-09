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

#ifndef OPENMS_ANALYSIS_OPENSWATH_SPECTRUMADDITON_H
#define OPENMS_ANALYSIS_OPENSWATH_SPECTRUMADDITON_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

namespace OpenMS
{
  /**
  @brief The SpectrumAddition adds together a list of spectra

  */
  class OPENMS_DLLAPI SpectrumAddition
  {

public:

    /// adds up a list of Spectra by resampling them and then addition of intensities
    static OpenSwath::SpectrumPtr addUpSpectra(std::vector<OpenSwath::SpectrumPtr> all_spectra, double sampling_rate, double filter_zeros)
    {
      if (all_spectra.size() == 0)
      {
        OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
        return sptr;
      }

      typedef MSSpectrum<Peak1D> SpectrumT;
      LinearResamplerAlign lresampler;
      // typedef std::vector<Peak1D> SpectrumT ;

      // find global min and max -> use as start/endpoints for resampling
      double min = all_spectra[0]->getMZArray()->data[0];
      double max = all_spectra[0]->getMZArray()->data.back();
      for (Size i = 0; i < all_spectra.size(); i++)
      {
        if (all_spectra[i]->getMZArray()->data[0] < min)
        {min = all_spectra[i]->getMZArray()->data[0]; }
        if (all_spectra[i]->getMZArray()->data.back() > max)
        {max = all_spectra[i]->getMZArray()->data.back(); }
      }

      // generate the resampled peaks at positions origin+i*spacing_
      int number_resampled_points = (max - min) / sampling_rate + 1;
      SpectrumT resampled_peak_container;
      resampled_peak_container.resize(number_resampled_points);
      SpectrumT::iterator it = resampled_peak_container.begin();
      for (int i = 0; i < number_resampled_points; ++i)
      {
        it->setMZ(min + i * sampling_rate);
        it->setIntensity(0);
        ++it;
      }

      // resample all spectra and add to master spectrum
      SpectrumT master_spectrum = resampled_peak_container;
      for (Size curr_sp = 0; curr_sp < all_spectra.size(); curr_sp++)
      {
        SpectrumT input_spectrum;
        SpectrumT output_spectrum = resampled_peak_container;

        // convert input spectrum to OpenMS, then resample
        OpenSwathDataAccessHelper::convertToOpenMSSpectrum(all_spectra[curr_sp], input_spectrum);
        lresampler.raster(input_spectrum.begin(), input_spectrum.end(), output_spectrum.begin(), output_spectrum.end());

        // add to master spectrum
        for (Size i = 0; i < output_spectrum.size(); ++i)
        {
          master_spectrum[i].setIntensity(master_spectrum[i].getIntensity() + output_spectrum[i].getIntensity());
        }
      }

      if (!filter_zeros)
      {
        OpenSwath::SpectrumPtr sptr = OpenSwathDataAccessHelper::convertToSpectrumPtr(master_spectrum);
        return sptr;
      }
      else
      {
        SpectrumT master_spectrum_filtered;
        for (Size i = 0; i < master_spectrum.size(); ++i)
        {
          if (master_spectrum[i].getIntensity() > 0)
          {
            master_spectrum_filtered.push_back(master_spectrum[i]);
          }
        }
        OpenSwath::SpectrumPtr sptr = OpenSwathDataAccessHelper::convertToSpectrumPtr(master_spectrum_filtered);
        return sptr;
      }
    }

  };
}

#endif
