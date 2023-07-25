// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

namespace OpenMS
{

  void sortSpectrumByMZ(OpenSwath::Spectrum& spec)
  {
    //sort index list
    std::vector<std::pair<double, Size> > sorted_indices;
    sorted_indices.reserve(spec.getMZArray()->data.size());
    auto mz_it = spec.getMZArray()->data.begin();
    for (Size i = 0; i < spec.getMZArray()->data.size(); ++i)
    {
      sorted_indices.emplace_back(*mz_it, i);
      ++mz_it;
    }
    std::stable_sort(sorted_indices.begin(), sorted_indices.end());

    // extract list of indices
    std::vector<Size> select_indices;
    select_indices.reserve(sorted_indices.size());
    for (const auto& sidx : sorted_indices)
    {
      select_indices.push_back(sidx.second);
    }

    for (auto& da : spec.getDataArrays() )
    {
      if (da->data.empty()) continue;
      OpenSwath::BinaryDataArrayPtr tmp(new OpenSwath::BinaryDataArray);
      tmp->description = da->description;
      tmp->data.reserve(select_indices.size());
      for (Size i = 0; i < select_indices.size(); ++i)
      {
        tmp->data.push_back( da->data[ select_indices[i] ] );
      }
      da = tmp;
    }

    OPENMS_POSTCONDITION( std::adjacent_find(spec.getMZArray()->data.begin(),
           spec.getMZArray()->data.end(), std::greater<double>()) == spec.getMZArray()->data.end(),
           "Postcondition violated: m/z vector needs to be sorted!" )
  }

  OpenSwath::SpectrumPtr SpectrumAddition::addUpSpectra(const SpectrumSequence& all_spectra, double sampling_rate, bool filter_zeros)
  {
    OPENMS_PRECONDITION(all_spectra.empty() || all_spectra[0]->getDataArrays().size() == 2, "Can only resample spectra with 2 data dimensions (no ion mobility spectra)")

    if (all_spectra.size() == 1) return all_spectra[0];
    if (all_spectra.empty())
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
    }
    // ensure first one is not empty
    if (all_spectra[0]->getMZArray()->data.empty() )
    {
      OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
      return sptr;
    }

    // find global min and max -> use as start/endpoints for resampling
    double min = all_spectra[0]->getMZArray()->data[0];
    double max = all_spectra[0]->getMZArray()->data.back();
    double min_spacing = max - min;
    for (Size i = 0; i < all_spectra.size(); i++)
    {
      if (all_spectra[i]->getMZArray()->data.empty() )
      {
        continue;
      }

      // estimate sampling rate
      for (Size k = 0; k < all_spectra[i]->getMZArray()->data.size() && sampling_rate < 0; k++)
      {
        if (k > 0)
        {
          if (min_spacing > all_spectra[i]->getMZArray()->data[k] - all_spectra[i]->getMZArray()->data[k-1] )
          {
            min_spacing = all_spectra[i]->getMZArray()->data[k] - all_spectra[i]->getMZArray()->data[k-1];
          }
        }
      }

      if (all_spectra[i]->getMZArray()->data[0] < min)
      {
        min = all_spectra[i]->getMZArray()->data[0];
      }
      if (all_spectra[i]->getMZArray()->data.back() > max)
      {
        max = all_spectra[i]->getMZArray()->data.back();
      }
    }

    // in case we are asked to estimate the resampling rate
    if (sampling_rate < 0) sampling_rate = min_spacing;

    // generate the resampled peaks at positions origin+i*spacing_
    int number_resampled_points = (max - min) / sampling_rate + 1;
    OpenSwath::SpectrumPtr resampled_peak_container(new OpenSwath::Spectrum);
    resampled_peak_container->getMZArray()->data.resize(number_resampled_points);
    resampled_peak_container->getIntensityArray()->data.resize(number_resampled_points);
    std::vector<double>::iterator it = resampled_peak_container->getMZArray()->data.begin();
    int cnt = 0;
    while (it != resampled_peak_container->getMZArray()->data.end())
    {
      *it = min + cnt * sampling_rate; // set mz (intensity is zero already)
      ++it;
      ++cnt;
    }

    LinearResamplerAlign lresampler;
    // resample all spectra and add to master spectrum
    for (Size curr_sp = 0; curr_sp < all_spectra.size(); curr_sp++)
    {
      lresampler.raster(all_spectra[curr_sp]->getMZArray()->data.begin(),
                        all_spectra[curr_sp]->getMZArray()->data.end(),
                        all_spectra[curr_sp]->getIntensityArray()->data.begin(),
                        all_spectra[curr_sp]->getIntensityArray()->data.end(),
                        resampled_peak_container->getMZArray()->data.begin(),
                        resampled_peak_container->getMZArray()->data.end(),
                        resampled_peak_container->getIntensityArray()->data.begin(),
                        resampled_peak_container->getIntensityArray()->data.end()
      );
    }

    if (!filter_zeros)
    {
      OPENMS_POSTCONDITION( std::adjacent_find(added_spec->getMZArray()->data.begin(),
           added_spec->getMZArray()->data.end(), std::greater<double>()) == added_spec->getMZArray()->data.end(),
           "Postcondition violated: m/z vector needs to be sorted!" )


      return resampled_peak_container;
    }
    else
    {
      OpenSwath::SpectrumPtr master_spectrum_filtered(new OpenSwath::Spectrum);
      for (Size i = 0; i < resampled_peak_container->getIntensityArray()->data.size(); ++i)
      {
        if (resampled_peak_container->getIntensityArray()->data[i] > 0)
        {
          master_spectrum_filtered->getIntensityArray()->data.push_back(resampled_peak_container->getIntensityArray()->data[i]);
          master_spectrum_filtered->getMZArray()->data.push_back(resampled_peak_container->getMZArray()->data[i]);
        }
      }

      OPENMS_POSTCONDITION( std::adjacent_find(added_spec->getMZArray()->data.begin(),
           added_spec->getMZArray()->data.end(), std::greater<double>()) == added_spec->getMZArray()->data.end(),
           "Postcondition violated: m/z vector needs to be sorted!" )


      return master_spectrum_filtered;
    }
  }


  OpenSwath::SpectrumPtr SpectrumAddition::addUpSpectra(const SpectrumSequence& all_spectra,
                                             const RangeMobility& im_range,
                                             double sampling_rate,
                                             bool filter_zeros)
  {
    OPENMS_PRECONDITION(! (im_range.isEmpty() || all_spectra[0].getFloatDataArrays().empty()), "Can only resample spectra with 2 data dimensions (no ion mobility spectra)")
    OpenSwath::SpectrumPtr added_spec(new OpenSwath::Spectrum);

    // If no spectra found return
    if (all_spectra.empty())
    {
      return added_spec;
    }

    if (im_range.isEmpty())
    {
      return addUpSpectra(all_spectra, sampling_rate, filter_zeros);
    }
    // since resampling is not supported on 3D data first filter by drift time (if possible) and then add
    // (!im_range.isEmpty())
    SpectrumSequence filteredSpectra;
    for (auto spec: all_spectra)
    {
      filteredSpectra.push_back(OpenSwath::ISpectrumAccess::filterByDrift(spec, im_range.getMin(), im_range.getMax()));
    }
    return addUpSpectra(filteredSpectra, sampling_rate, filter_zeros);
  }

  OpenSwath::SpectrumPtr SpectrumAddition::concatenateSpectra(const SpectrumSequence& all_spectra)

  {
    OpenSwath::SpectrumPtr added_spec(new OpenSwath::Spectrum);
    // Ensure that we have the same number of data arrays as in the input spectrum
    // copying the extra data arrays descriptions onto the added spectra
    if (!all_spectra.empty() && all_spectra[0]->getDataArrays().size() > 2)
    {
      for (Size k = 2; k < all_spectra[0]->getDataArrays().size(); k++)
      {
        OpenSwath::BinaryDataArrayPtr tmp (new OpenSwath::BinaryDataArray());
        tmp->description = all_spectra[0]->getDataArrays()[k]->description;
        added_spec->getDataArrays().push_back(tmp);
      }
    }

    // Simply concatenate all spectra together and sort in the end
    for (const auto& s : all_spectra)
    {
      for (Size k = 0; k < s->getDataArrays().size(); k++)
      {
        auto& v1 = added_spec->getDataArrays()[k]->data;
        auto& v2 = s->getDataArrays()[k]->data;

        v1.reserve( v1.size() + v2.size() );
        v1.insert( v1.end(), v2.begin(), v2.end() );
      }
    }
    sortSpectrumByMZ(*added_spec);
    return added_spec;
  }

  OpenMS::MSSpectrum SpectrumAddition::addUpSpectra(const std::vector<MSSpectrum>& all_spectra, double sampling_rate, bool filter_zeros)
  {
    OPENMS_PRECONDITION(all_spectra.empty() || all_spectra[0].getFloatDataArrays().empty(), "Can only resample spectra with 2 data dimensions (no ion mobility spectra)")

    if (all_spectra.size() == 1) return all_spectra[0];
    if (all_spectra.empty()) return MSSpectrum();
    // ensure first one is not empty
    if (all_spectra[0].empty() ) return MSSpectrum();

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

      // estimate sampling rate
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

    // in case we are asked to estimate the resampling rate
    if (sampling_rate < 0) sampling_rate = min_spacing;

    // generate the resampled peaks at positions origin+i*spacing_
    int number_resampled_points = (max - min) / sampling_rate + 1;
    MSSpectrum resampled_peak_container;
    resampled_peak_container.resize(number_resampled_points);
    MSSpectrum::iterator it = resampled_peak_container.begin();
    for (int i = 0; i < number_resampled_points; ++i)
    {
      it->setMZ(min + i * sampling_rate);
      it->setIntensity(0);
      ++it;
    }

    // resample all spectra and add to master spectrum
    LinearResamplerAlign lresampler;
    MSSpectrum master_spectrum = resampled_peak_container;
    for (Size curr_sp = 0; curr_sp < all_spectra.size(); curr_sp++)
    {
      MSSpectrum input_spectrum;
      MSSpectrum output_spectrum = resampled_peak_container;

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
      MSSpectrum master_spectrum_filtered;
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
