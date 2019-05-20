// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

// #define DEBUG

using namespace std;
using namespace boost::math;

namespace OpenMS
{

  MultiplexFilteringCentroided::MultiplexFilteringCentroided(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_centroided, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type)
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
        for (MSSpectrum::ConstIterator it_mz = it_rt.begin(); it_mz != it_rt.end(); ++it_mz)
        {
          double mz = it_mz->getMZ();
          MultiplexFilteredPeak peak(mz, rt, exp_centroided_mapping_[idx_rt][it_mz - it_rt.begin()], idx_rt);
          
          if (!(filterPeakPositions_(it_mz, exp_centroided_white_.begin(), it_rt_band_begin, it_rt_band_end, pattern, peak)))
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

          result.addPeak(peak);
          blacklistPeak_(peak, pattern_idx);
        }
      }
      
#ifdef DEBUG
      // write filtered peaks to debug output
      std::stringstream debug_out;
      debug_out << "filter_result_" << pattern_idx << ".consensusXML";
      result.writeDebugOutput(exp_centroided_, debug_out.str());
#endif
      
      // add results of this pattern to list
      filter_results.push_back(result);
    }
    
#ifdef DEBUG
    // clock for monitoring run time performance
    OPENMS_LOG_INFO << "\nThe filtering step of the algorithm took " << (float)(clock()-start)/CLOCKS_PER_SEC << " seconds.\n\n";
#endif

    endProgress();
    
    return filter_results;
  }
  
}
