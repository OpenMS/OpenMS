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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
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

  MultiplexFilteringCentroided::MultiplexFilteringCentroided(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type)
  {
  }

  vector<MultiplexFilteredMSExperiment> MultiplexFilteringCentroided::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_picked_.size(), "filtering LC-MS data");
    
    // list of filter results for each peak pattern
    vector<MultiplexFilteredMSExperiment> filter_results;
    
    unsigned int start = clock();
    
    // loop over all patterns
    for (unsigned pattern_idx = 0; pattern_idx < patterns_.size(); ++pattern_idx)
    {
      // current pattern
      MultiplexIsotopicPeakPattern pattern = patterns_[pattern_idx];
      
      // data structure storing peaks which pass all filters for this pattern
      MultiplexFilteredMSExperiment result;
  
      // construct new white experiment
      White2Original exp_picked_mapping;
      MSExperiment exp_picked_white = getWhiteMSExperiment_(exp_picked_mapping);
  
      // filter (white) experiment
      // loop over spectra
      for (MSExperiment::ConstIterator it_rt = exp_picked_white.begin(); it_rt < exp_picked_white.end(); ++it_rt)
      {
        // skip empty spectra
        if (it_rt->empty())
        {
          continue;
        }

        setProgress(++progress);

        double rt = it_rt->getRT();
        MSExperiment::ConstIterator it_rt_band_begin = exp_picked_white.RTBegin(rt - rt_band_/2);
        MSExperiment::ConstIterator it_rt_band_end = exp_picked_white.RTEnd(rt + rt_band_/2);
        
        // debug output variables
        /*size_t debug_pattern_idx = 8;
        size_t debug_rt_idx = 31;
        size_t debug_mz_idx = 4;*/
        
        // debug output
        /*if (pattern_idx == debug_pattern_idx)
        {
         std::cout << "RT = " << rt << "    RT idx = " << (it_rt - exp_picked_white.begin()) << "\n";
        }*/
        
        // loop over m/z
        for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
        {
          double mz = it_mz->getMZ();
          MultiplexFilteredPeak peak(mz, rt, exp_picked_mapping[it_rt - exp_picked_white.begin()][it_mz - it_rt->begin()], it_rt - exp_picked_white.begin());
          
          // debug output
          //bool debug_now = (pattern_idx == debug_pattern_idx) && (it_rt - exp_picked_white.begin() == debug_rt_idx) && (it_mz - it_rt->begin() == debug_mz_idx);
          /*if ((pattern_idx == debug_pattern_idx) && (it_rt - exp_picked_white.begin() == debug_rt_idx))
          {
            std::cout << "RT = " << rt << "    m/z = " << mz << "    m/z idx = " << (it_mz - it_rt->begin()) << "    int = " << it_mz->getIntensity() << "\n";
          }*/
          /*if (debug_now)
          {
            std::cout << "\n";
            std::cout << "RT = " << rt << "    m/z = " << mz << "    int = " << it_mz->getIntensity() << "\n";
          }*/

          if (!(filterPeakPositions_(it_mz, exp_picked_mapping, exp_picked_white.begin(), it_rt_band_begin, it_rt_band_end, pattern, peak)))
          {
            continue;
          }
          
          // debug output
          /*if (debug_now)
          {
            std::cout << "Passed Peak Position Filter.\n";
          }*/

          if (!(filterAveragineModel_(pattern, peak)))
          {
            continue;
          }
          
          // debug output
          /*if (debug_now)
          {
            std::cout << "Passed Averagine Filter.\n";
          }*/

          if (!(filterPeptideCorrelation_(pattern, peak)))
          {
            continue;
          }
                    
          // debug output
          /*if (debug_now)
          {
            std::cout << "Passed peptide similarity filter.\n\n\n";
          }*/

          result.addPeak(peak);
          //blacklistPeak_(peak);
          blacklistPeak2_(peak, pattern_idx);
        }
      }
      
      // write filtered peaks to debug output
      std::stringstream debug_out;
      debug_out << "filter_result_" << pattern_idx << ".consensusXML";
      result.writeDebugOutput(exp_picked_, debug_out.str());
      
      // add results of this pattern to list
      filter_results.push_back(result);
      
      ungreyBlacklist_();
    }
            
    std::cout << "\nFiltering took me " << (float)(clock()-start)/CLOCKS_PER_SEC << " seconds.\n\n";

    endProgress();
    
    return filter_results;
  }
  
}
