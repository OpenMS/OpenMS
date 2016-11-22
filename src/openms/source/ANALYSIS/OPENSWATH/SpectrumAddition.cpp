// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

namespace OpenMS
{

  OpenSwath::SpectrumPtr SpectrumAddition::addUpSpectra(std::vector<OpenSwath::SpectrumPtr> all_spectra,
      double sampling_rate, bool filter_zeros)
  {
    if (all_spectra.empty())
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
    }

    typedef MSSpectrum<Peak1D> SpectrumT;
    LinearResamplerAlign lresampler;

    if (all_spectra[0]->getMZArray()->data.empty() )
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
    }

    // find global min and max -> use as start/endpoints for resampling
    double min = all_spectra[0]->getMZArray()->data[0];
    double max = all_spectra[0]->getMZArray()->data.back();
    bool all_empty = true;
    for (Size i = 0; i < all_spectra.size(); i++)
    {
      if (all_spectra[i]->getMZArray()->data.empty() )
      {
        continue;
      }
      all_empty = false;

      if (all_spectra[i]->getMZArray()->data[0] < min)
      {
        min = all_spectra[i]->getMZArray()->data[0];
      }
      if (all_spectra[i]->getMZArray()->data.back() > max)
      {
        max = all_spectra[i]->getMZArray()->data.back();
      }
    }

    if (all_empty)
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
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


  OpenMS::MSSpectrum<> SpectrumAddition::addUpSpectra(std::vector<OpenMS::MSSpectrum<> > all_spectra, double sampling_rate, bool filter_zeros)
  {
    if (all_spectra.size() == 1)
    {
      return all_spectra[0];
    }

    bool all_empty = true;
    for (Size i = 0; i < all_spectra.size(); i++)
    {
      if (all_spectra[i].empty())
      {
        continue;
      }
      all_empty = false;
    }
    if (all_spectra.empty() || all_empty)
    {
      return MSSpectrum<>();
    }
    if (all_spectra[0].empty() )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "First spectrum cannot be empty");
    }

    // find global min and max -> use as start/endpoints for resampling
    double min = all_spectra[0][0].getMZ();
    double max = all_spectra[0][all_spectra[0].size()-1].getMZ();
    double min_spacing = max - min;
    for (Size i = 0; i < all_spectra.size(); i++)
    {
      if (all_spectra[i].empty())
      {
        continue;
      }

      for (Size k = 0; k < all_spectra[i].size() && sampling_rate < 0; k++)
      {
        if (k > 0)
        {
          if (min_spacing > all_spectra[i][k].getMZ() - all_spectra[i][k-1].getMZ() )
          {
            min_spacing = all_spectra[i][k].getMZ() - all_spectra[i][k-1].getMZ();
          }
        }

        if (all_spectra[i][k].getMZ() < min) min = all_spectra[i][k].getMZ();
        if (all_spectra[i][k].getMZ() > max) max = all_spectra[i][k].getMZ();
      }

      if (all_spectra[i][0].getMZ() < min) min = all_spectra[i][0].getMZ();
      if (all_spectra[i][ all_spectra[i].size() -1].getMZ() > max) max = all_spectra[i][ all_spectra[i].size() -1].getMZ();
    }

    if (all_empty)
    {
      return MSSpectrum<>();
    }

    typedef MSSpectrum<Peak1D> SpectrumT;
    LinearResamplerAlign lresampler;

    // in case we are asked to estimate the resampling rate
    if (sampling_rate < 0) sampling_rate = min_spacing;

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

      lresampler.raster(all_spectra[curr_sp].begin(), all_spectra[curr_sp].end(), output_spectrum.begin(), output_spectrum.end());

      // add to master spectrum
      for (Size i = 0; i < output_spectrum.size(); ++i)
      {
        master_spectrum[i].setIntensity(master_spectrum[i].getIntensity() + output_spectrum[i].getIntensity());
      }
    }

    if (!filter_zeros)
    {
      return master_spectrum;
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
      return master_spectrum_filtered;
    }
  }

}

