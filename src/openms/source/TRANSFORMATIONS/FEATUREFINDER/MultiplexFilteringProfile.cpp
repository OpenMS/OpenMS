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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <QDir>

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFilteringProfile::MultiplexFilteringProfile(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, const std::vector<MultiplexIsotopicPeakPattern> patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type), exp_profile_(exp_profile), boundaries_(boundaries)
  {
    
    if (exp_profile_.size() != exp_picked_.size())
    {
      stringstream stream;
      stream << "Profile and centroided data do not contain same number of spectra. (";
      stream << exp_profile_.size();
      stream << "!=";
      stream << exp_picked_.size();
      stream << ")";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
    }

    if (exp_picked_.size() != boundaries_.size())
    {
      stringstream stream;
      stream << "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra. (";
      stream << exp_picked_.size();
      stream << "!=";
      stream << boundaries_.size();
      stream << ")";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
    }

  }

  vector<MultiplexFilteredMSExperiment> MultiplexFilteringProfile::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_profile_.size(), "filtering LC-MS data");

    // list of filter results for each peak pattern
    std::vector<MultiplexFilteredMSExperiment> filter_results;
    
    std::cout << "\nStart filtering.\n\n";
      
    unsigned int start = clock();
    
    // loop over all patterns
    for (unsigned pattern_idx = 0; pattern_idx < patterns_.size(); ++pattern_idx)
    {
      std::cout << "\npattern = " << pattern_idx << "\n";
      
      // current pattern
      MultiplexIsotopicPeakPattern pattern = patterns_[pattern_idx];
      
      // data structure storing peaks which pass all filters
      MultiplexFilteredMSExperiment result;
  
      // loop over spectra
      MSExperiment::ConstIterator it_rt_picked;
      MSExperiment::ConstIterator it_rt_profile;
      std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
      for (it_rt_picked = exp_picked_.begin(), it_rt_profile = exp_profile_.begin(), it_rt_boundaries = boundaries_.begin();
           it_rt_picked < exp_picked_.end() && it_rt_profile < exp_profile_.end() && it_rt_boundaries < boundaries_.end();
           ++it_rt_picked, ++it_rt_profile, ++it_rt_boundaries)
      {
        // skip empty spectra
        if ((*it_rt_profile).size() == 0 || (*it_rt_picked).size() == 0 || (*it_rt_boundaries).size() == 0)
        {
          continue;
        }
        
        if ((*it_rt_picked).size() != (*it_rt_boundaries).size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Number of peaks and number of peak boundaries differ.");
        }

        setProgress(++progress);
        
        double rt = it_rt_picked->getRT();
        double rt_min = rt - rt_band_/2;
        double rt_max = rt + rt_band_/2;
        
        std::cout << "    RT = " << rt << "    RT band (min) = " << rt_min << "    RT band (max) = " << rt_max << "\n";
        
        // loop over mz
        for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt_picked->begin(); it_mz < it_rt_picked->end(); ++it_mz)
        {
          double mz = it_mz->getMZ();
          std::cout << "        mz = " << mz << "\n";
          
          // loop over peaks in pattern
          /*for (size_t mz_position = 1; mz_position < pattern.getMZShiftCount(); ++mz_position)
          {
            double mz_shift = pattern.getMZShiftAt(mz_position);
          
            // loop over RT band
            int count = 0;
            double intensity = 0;
            for (MSExperiment::ConstIterator it_rt_band = exp_picked_.RTBegin(rt_min); it_rt_band < exp_picked_.RTEnd(rt_max); ++it_rt_band)
            {
              int i = it_rt_band->findNearest(mz + mz_shift, mz_tol);
              
              if (i == -1)
              {
                continue;
              }
              
              ++count;
              //intensity += it_rt[i].getMZ();
              
              std::cout << "        RT = " << it_rt_band->getRT() << "        i = " << i << "\n";
            }
            
            intensity = intensity/count;
            
          }*/
        }
      }
     
      // add results of this pattern to list
      filter_results.push_back(result);
    }
        
    std::cout << "\nThat took me " << (float)(clock()-start)/CLOCKS_PER_SEC << " seconds.\n";
    std::cout << "\nFinished filtering.\n\n";

    endProgress();

    return filter_results;
  }

  int MultiplexFilteringProfile::nonLocalIntensityFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& mz_shifts_actual, const vector<int>& mz_shifts_actual_indices, SplineSpectrum::Navigator nav, std::vector<double>& intensities_actual, int peaks_found_in_all_peptides, double mz) const
  {
    // calculate intensities
    for (int i = 0; i < (int) mz_shifts_actual_indices.size(); ++i)
    {
      if (mz_shifts_actual_indices[i] != -1)
      {
        intensities_actual.push_back(nav.eval(mz + mz_shifts_actual[i]));
      }
      else
      {
        intensities_actual.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }

    // peaks_found_in_all_peptides is the number of peaks in this region. At this particular m/z
    // some of the intensities might be below the cutoff. => return isotope
    for (int isotope = 0; isotope < peaks_found_in_all_peptides; ++isotope)
    {
      bool seen_in_all_peptides = true;
      for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
      {
        if (boost::math::isnan(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1]))
        {
          // peak not found
          seen_in_all_peptides = false;
          break;
        }
        else if (intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1] < intensity_cutoff_)
        {
          // below intensity threshold
          seen_in_all_peptides = false;
          break;
        }
      }
      if (!seen_in_all_peptides)
      {
        return isotope;
      }
    }

    return peaks_found_in_all_peptides;
  }

  int MultiplexFilteringProfile::findNearest_(int spectrum_index, double mz, double scaling) const
  {
    MSExperiment::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries = boundaries_.begin() + spectrum_index;

    MSSpectrum<Peak1D>::ConstIterator it_mz;
    vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundaries;
    for (it_mz = it_rt->begin(), it_mz_boundaries = it_rt_boundaries->begin();
         it_mz < it_rt->end() && it_mz_boundaries < it_rt_boundaries->end();
         ++it_mz, ++it_mz_boundaries)
    {
      if (mz >= scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ() &&
          mz <= scaling * (*it_mz_boundaries).mz_max + (1 - scaling) * it_mz->getMZ())
      {
        return it_mz - it_rt->begin();
      }
      if (mz < scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ())
      {
        return -1;
      }
    }

    return -1;
  }

}
