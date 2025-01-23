// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

// #define DEBUG

using namespace std;

namespace OpenMS
{

  MultiplexFilteringCentroided::MultiplexFilteringCentroided(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_centroided, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, std::move(averagine_type))
  {
  }

  vector<MultiplexFilteredMSExperiment> MultiplexFilteringCentroided::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_centroided_.size(), "filtering LC-MS data");
    
    // list of filter results for each peak pattern
    vector<MultiplexFilteredMSExperiment> filter_results;

#ifdef DEBUG
    // clock for monitoring run time performance
    unsigned int start = clock();
#endif

    // loop over all patterns
    for (unsigned pattern_idx = 0; pattern_idx < patterns_.size(); ++pattern_idx)
    {
      // current pattern
      MultiplexIsotopicPeakPattern pattern = patterns_[pattern_idx];
      
      // data structure storing peaks which pass all filters for this pattern
      MultiplexFilteredMSExperiment result;
  
      // update white experiment
      updateWhiteMSExperiment_();
  
      // filter (white) experiment
      // loop over spectra
      for (const auto &it_rt : exp_centroided_white_)
      {
        // skip empty spectra
        if (it_rt.empty())
        {
          continue;
        }

        setProgress(++progress);

        double rt = it_rt.getRT();
        size_t idx_rt = &it_rt - &exp_centroided_white_[0];
        
        MSExperiment::ConstIterator it_rt_band_begin = exp_centroided_white_.RTBegin(rt - rt_band_/2);
        MSExperiment::ConstIterator it_rt_band_end = exp_centroided_white_.RTEnd(rt + rt_band_/2);
        
        // loop over m/z
        //for (MSSpectrum::ConstIterator it_mz = it_rt.begin(); it_mz != it_rt.end(); ++it_mz)
        #pragma omp parallel for
        for (SignedSize s = 0; s < (SignedSize) it_rt.size(); s++)
        {
          auto& it_mz = it_rt[s];
          double mz = it_mz.getMZ();
          MultiplexFilteredPeak peak(mz, rt, exp_centroided_mapping_[idx_rt][s], idx_rt);
          
          if (!(filterPeakPositions_(mz, exp_centroided_white_.begin(), it_rt_band_begin, it_rt_band_end, pattern, peak)))
          {
            continue;
          }
          
          if (!(filterAveragineModel_(pattern, peak)))
          {
            continue;
          }

          if (!(filterPeptideCorrelation_(pattern, peak)))
          {
            continue;
          }
          
          /**
           * All filters passed.
           */

          #pragma omp critical
          {
            result.addPeak(peak);
            blacklistPeak_(peak, pattern_idx);
          };
        }
      }
      
      // add results of this pattern to list
      filter_results.push_back(result);
    }
    
    endProgress();
    
    return filter_results;
  }
  
}
