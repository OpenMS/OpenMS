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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/STATISTICS/LinearRegressionWithoutIntercept.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>

// #define DEBUG

using namespace std;

namespace OpenMS
{
  FeatureFinderMultiplexAlgorithm::FeatureFinderMultiplexAlgorithm() :
    DefaultParamHandler("FeatureFinderMultiplexAlgorithm")
  {
    // parameter section: algorithm
    defaults_.setValue("algorithm:labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. If the sample is unlabelled (i.e. you want to detect only single peptide features) please leave this parameter empty. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
    defaults_.setValue("algorithm:charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
    defaults_.setValue("algorithm:isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", ListUtils::create<String>("advanced"));
    defaults_.setValue("algorithm:rt_typical", 40.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
    defaults_.setMinFloat("algorithm:rt_typical", 0.0);
    defaults_.setValue("algorithm:rt_band", 10.0, "RT band which is taken into considerations when filtering.TODO docu");
    defaults_.setMinFloat("algorithm:rt_band", 0.0);
    defaults_.setValue("algorithm:rt_min", 2.0, "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)");
    defaults_.setMinFloat("algorithm:rt_min", 0.0);
    defaults_.setValue("algorithm:mz_tolerance", 6.0, "m/z tolerance for search of peak patterns.");
    defaults_.setMinFloat("algorithm:mz_tolerance", 0.0);
    defaults_.setValue("algorithm:mz_unit", "ppm", "Unit of the 'mz_tolerance' parameter.");
    defaults_.setValidStrings("algorithm:mz_unit", ListUtils::create<String>("Da,ppm"));
    defaults_.setValue("algorithm:intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks.");
    defaults_.setMinFloat("algorithm:intensity_cutoff", 0.0);
    defaults_.setValue("algorithm:peptide_similarity", 0.5, "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.");
    defaults_.setMinFloat("algorithm:peptide_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:peptide_similarity", 1.0);
    defaults_.setValue("algorithm:averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
    defaults_.setMinFloat("algorithm:averagine_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity", 1.0);
    defaults_.setValue("algorithm:averagine_similarity_scaling", 0.95, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simulataneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("algorithm:averagine_similarity_scaling", 0.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity_scaling", 1.0);
    defaults_.setValue("algorithm:missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
    defaults_.setMinInt("algorithm:missed_cleavages", 0);
    defaults_.setValue("algorithm:spectrum_type", "automatic", "Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("algorithm:spectrum_type", ListUtils::create<String>("profile,centroid,automatic"));
    defaults_.setValue("algorithm:averagine_type","peptide","The type of averagine to use, currently RNA, DNA or peptide", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("algorithm:averagine_type", ListUtils::create<String>("peptide,RNA,DNA"));
    defaults_.setValue("algorithm:knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("algorithm:knock_out", ListUtils::create<String>("true,false"));

    defaults_.setSectionDescription("algorithm", "algorithmic parameters");
    
    // parameter section: labels
    defaults_.setValue("labels:Arg6", 12.08, "description", ListUtils::create<String>("advanced"));
    
    defaults_.setSectionDescription("labels", "mass shifts for all possible labels");
    
    MultiplexDeltaMassesGenerator generator;
    Param p = generator.getParameters();
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      String label_name = "labels:";
      label_name += it->name;
      
      defaults_.setValue(label_name, it->value, it->description, ListUtils::create<String>("advanced"));
      defaults_.setMinFloat(label_name, 0.0);
      
      label_mass_shift_.insert(make_pair(it->name, it->value));
    }
    
    // parameter section: algorithm, get selected charge range
    String charge_string = defaults_.getValue("algorithm:charge");
    charge_min_ = charge_string.prefix(':').toInt();
    charge_max_ = charge_string.suffix(':').toInt();
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }
    
    // parameter section: algorithm, get isotopes per peptide range
    String isotopes_per_peptide_string = defaults_.getValue("algorithm:isotopes_per_peptide");
    isotopes_per_peptide_min_ = isotopes_per_peptide_string.prefix(':').toInt();
    isotopes_per_peptide_max_ = isotopes_per_peptide_string.suffix(':').toInt();
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }
  }
  
  /**
   * @brief order of charge states
   *
   * 2+ 3+ 4+ 1+ 5+ 6+ ...
   *
   * Order charge states by the likelihood of their occurrence, i.e. we search for the most likely charge states first.
   */
  static size_t orderCharge(int charge)
  {
    if ((1 < charge) && (charge < 5))
    {
      return (charge - 1);
    }
    else if (charge == 1)
    {
      return 4;
    }
    else
    {
      return charge;
    }
  }

  /**
   * @brief comparator of peak patterns
   *
   * The comperator determines in which order the peak patterns are searched for.
   * First we check the number of mass shifts (triplets before doublets before singlets).
   * Then we check the first mass shift (for example 6 Da before 12 Da i.e. misscleavage).
   * Finally we check for charges (2+ before 1+, most likely first).
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool lessPattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      // The first mass shift is by definition always zero.
      if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
      {
        if (pattern1.getMassShiftAt(1) == pattern2.getMassShiftAt(1))
        {
          // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
          return orderCharge(pattern1.getCharge()) < orderCharge(pattern2.getCharge());
        }
        else
        {
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
      }
      else
      {
        // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
        return orderCharge(pattern1.getCharge()) < orderCharge(pattern2.getCharge());
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  std::vector<MultiplexIsotopicPeakPattern> FeatureFinderMultiplexAlgorithm::generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, const std::vector<MultiplexDeltaMasses>& mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;
    
    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(), list.end(), lessPattern);
    
#ifdef DEBUG
    // debug output
    for (int i = 0; i < list.size(); ++i)
    {
      std::cout << "charge = " << list[i].getCharge() << "+    shift = " << list[i].getMassShiftAt(1) << " Da\n";
    }
#endif
    
    return list;
  }
  
  std::vector<double> FeatureFinderMultiplexAlgorithm::determinePeptideIntensitiesCentroided_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteCentroided >& satellites)
  {
    // determine RT shift between the peptides
    // i.e. first determine the RT centre of mass for each peptide
    std::vector<double> rt_peptide;
    std::vector<double> intensity_peptide;
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // coordinates of the peptide feature
      // RT is the intensity-average of all satellites peaks of all (!) mass traces
      double rt(0);
      double intensity_sum(0);
      
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope i.e. mass trace
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope;
        satellites_isotope = satellites.equal_range(idx);
        
        // loop over satellites for this isotope i.e. mass trace
        for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
        {
          // find indices of the peak
          size_t rt_idx = (satellite_it->second).getRTidx();
          size_t mz_idx = (satellite_it->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
          std::advance(it_rt, rt_idx);
          MSSpectrum::ConstIterator it_mz = it_rt->begin();
          std::advance(it_mz, mz_idx);
          
          rt += it_rt->getRT() * it_mz->getIntensity();
          intensity_sum += it_mz->getIntensity();
        }
      }
      
      rt /= intensity_sum;
      rt_peptide.push_back(rt);
      intensity_peptide.push_back(intensity_sum);
    }
    
    // determine the fold changes between the lightest peptide and the remaining ones
    // TODO Replace the lightest peptide by the highest intensity peptide in the multiplet
    std::vector<double> ratio_peptide;
    ratio_peptide.push_back(1.0);
    // loop over other peptides
    for (size_t peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // fill the vectors from which the Pearson correlation for the ratio i.e. fold change will be calculated
      std::vector<double> intensities_light;
      std::vector<double> intensities_other;
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope in the light peptide
        std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope_1;
        satellites_isotope_1 = satellites.equal_range(isotope);
        
        // find satellites for this isotope in the second peptide
        std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope_2;
        satellites_isotope_2 = satellites.equal_range(peptide * isotopes_per_peptide_max_ + isotope);
        
        // loop over satellites for this isotope in the light peptide
        for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_1 = satellites_isotope_1.first; satellite_it_1 != satellites_isotope_1.second; ++satellite_it_1)
        {
          // find indices of the peak
          size_t rt_idx_1 = (satellite_it_1->second).getRTidx();
          size_t mz_idx_1 = (satellite_it_1->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt_1 = exp_centroid_.begin();
          std::advance(it_rt_1, rt_idx_1);
          MSSpectrum::ConstIterator it_mz_1 = it_rt_1->begin();
          std::advance(it_mz_1, mz_idx_1);
          
          // find corresponding spectra
          double rt_1 = it_rt_1->getRT();
          double rt_2_target = rt_1 + rt_peptide[peptide] - rt_peptide[0];
          
          MSExperiment::ConstIterator it_rt_2_tmp = exp_centroid_.RTBegin(rt_2_target);
          double rt_2 = it_rt_2_tmp->getRT();
          // The previous spectrum might be a better match for the target RT <rt_2_target>.
          if (it_rt_2_tmp != exp_centroid_.begin())
          {
            if (std::abs((it_rt_2_tmp - 1)->getRT() - rt_2_target) < std::abs(rt_2 - rt_2_target))
            {
              --it_rt_2_tmp;
              rt_2 = it_rt_2_tmp->getRT();
            }
          }
          
          // loop over satellites for this isotope in the second peptide
          double rt_earlier = -1;
          double intensity_earlier = -1;
          double rt_later = -1;
          double intensity_later = -1;
          for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_2 = satellites_isotope_2.first; satellite_it_2 != satellites_isotope_2.second; ++satellite_it_2)
          {
            // find indices of the peak
            size_t rt_idx_2 = (satellite_it_2->second).getRTidx();
            size_t mz_idx_2 = (satellite_it_2->second).getMZidx();
            
            // find peak itself
            MSExperiment::ConstIterator it_rt_2 = exp_centroid_.begin();
            std::advance(it_rt_2, rt_idx_2);
            MSSpectrum::ConstIterator it_mz_2 = it_rt_2->begin();
            std::advance(it_mz_2, mz_idx_2);
            
            if (it_rt_2->getRT() <= rt_2 && (std::abs(it_rt_2->getRT() - rt_2) < std::abs(rt_earlier - rt_2)))
            {
              rt_earlier = it_rt_2->getRT();
              intensity_earlier = it_mz_2->getIntensity();
            }
            
            if (it_rt_2->getRT() >= rt_2 && (std::abs(it_rt_2->getRT() - rt_2) < std::abs(rt_later - rt_2)))
            {
              rt_later = it_rt_2->getRT();
              intensity_later = it_mz_2->getIntensity();
            }
            
          }
          
          // Our target lies on or between two satellites of the 'other' peptide.
          if ((rt_earlier > 0) && (rt_later > 0))
          {
            // linearly interpolated intensity
            double intensity_other;
            if ((rt_2 == rt_earlier) || (rt_later == rt_earlier))
            {
              intensity_other = intensity_earlier;
            }
            else
            {
              intensity_other = intensity_earlier + (intensity_later - intensity_earlier)*(rt_2 - rt_earlier)/(rt_later - rt_earlier);
            }
            
            intensities_light.push_back(it_mz_1->getIntensity());
            intensities_other.push_back(intensity_other);
          }
        }
        
      }
      
      // If less than three matches are found, we cannot reliably calculate the intensity ratio (aka slope) and report the uncorrected intensities.
      if ((intensities_light.size() < 3) || (intensities_other.size() < 3))
      {
        return intensity_peptide;
      }
      
      // determine ratios through linear regression of all corresponding intensities
      Math::LinearRegressionWithoutIntercept linreg;
      linreg.addData(intensities_light, intensities_other);
      double slope = linreg.getSlope();
      
      ratio_peptide.push_back(slope);
    }
    
    // correct peptide intensities
    // The peptide ratios are calculated as linear regression of (spline-interpolated) profile intensities, @see linreg
    // The individual peptide intensities are the sum of the same profile intensities. But the quotient of these peptide intensities
    // is not necessarily the same as the independently calculated ratio from the linear regression. Since the peptide ratio
    // from linear regression is the more accurate one, we correct the two peptide intensities by projecting them onto the ratio.
    // In the end, both peptide ratio from linear regression and the quotient of the peptide intensities are identical.
    std::vector<double> intensity_peptide_corrected;
    if (intensity_peptide.size() == 2)
    {
      double intensity1 = (intensity_peptide[0] + ratio_peptide[1] * intensity_peptide[1]) / (1 + ratio_peptide[1] * ratio_peptide[1]);
      double intensity2 = ratio_peptide[1] * intensity1;
      intensity_peptide_corrected.push_back(intensity1);
      intensity_peptide_corrected.push_back(intensity2);
    }
    else if (intensity_peptide.size() > 2)
    {
      // Now with n instead of two peptide intensities, one needs to project the peptide intensities onto the hyperplane defined
      // by the set of all peptide ratios (TODO). Instead, it is simpler to keep the lightest peptide intensity fixed, and correct
      // only the remaining ones. The correct peptide ratio (from linear regression) is reported on both cases.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
      for (unsigned i = 1; i < intensity_peptide.size(); ++i)
      {
        intensity_peptide_corrected.push_back(ratio_peptide[i] * intensity_peptide[0]);
      }
    }
    else
    {
      // For simple feature detection (singlets) the intensities remain unchanged.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
    }
    
    return intensity_peptide_corrected;
  }

  std::vector<double> FeatureFinderMultiplexAlgorithm::determinePeptideIntensitiesProfile_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites)
  {
    // determine RT shift between the peptides
    // i.e. first determine the RT centre of mass for each peptide
    std::vector<double> rt_peptide;
    std::vector<double> intensity_peptide;
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // coordinates of the peptide feature
      // RT is the intensity-average of all satellites peaks of all (!) mass traces
      double rt(0);
      double intensity_sum(0);
      
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope i.e. mass trace
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope;
        satellites_isotope = satellites.equal_range(idx);
        
        // loop over satellites for this isotope i.e. mass trace
        for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
        {
          rt += (satellite_it->second).getRT() * (satellite_it->second).getIntensity();
          intensity_sum += (satellite_it->second).getIntensity();
        }
      }
      
      rt /= intensity_sum;
      rt_peptide.push_back(rt);
      intensity_peptide.push_back(intensity_sum);
    }
    
    // determine the fold changes between the lightest peptide and the remaining ones
    // TODO Replace the lightest peptide by the highest intensity peptide in the multiplet
    std::vector<double> ratio_peptide;
    ratio_peptide.push_back(1.0);
    // loop over other peptides
    for (size_t peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // fill the vectors from which the Pearson correlation for the ratio i.e. fold change will be calculated
      std::vector<double> intensities_light;
      std::vector<double> intensities_other;
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope in the light peptide
        std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope_1;
        satellites_isotope_1 = satellites.equal_range(isotope);
        
        // find satellites for this isotope in the second peptide
        std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope_2;
        satellites_isotope_2 = satellites.equal_range(peptide * isotopes_per_peptide_max_ + isotope);
        
        // loop over satellites for this isotope in the light peptide
        for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it_1 = satellites_isotope_1.first; satellite_it_1 != satellites_isotope_1.second; ++satellite_it_1)
        {
          // find corresponding spectra
          double rt_1 = (satellite_it_1->second).getRT();
          double rt_2_target = rt_1 + rt_peptide[peptide] - rt_peptide[0];
          
          MSExperiment::ConstIterator it_rt_2 = exp_centroid_.RTBegin(rt_2_target);
          double rt_2 = it_rt_2->getRT();
          // The previous spectrum might be a better match for the target RT <rt_2_target>.
          if (it_rt_2 != exp_centroid_.begin())
          {
            if (std::abs((it_rt_2 - 1)->getRT() - rt_2_target) < std::abs(rt_2 - rt_2_target))
            {
              --it_rt_2;
              rt_2 = it_rt_2->getRT();
            }
          }
          
          // loop over satellites for this isotope in the second peptide
          double rt_earlier = -1;
          double intensity_earlier = -1;
          double rt_later = -1;
          double intensity_later = -1;
          for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it_2 = satellites_isotope_2.first; satellite_it_2 != satellites_isotope_2.second; ++satellite_it_2)
          {
            if ((satellite_it_2->second).getRT() <= rt_2 && (std::abs((satellite_it_2->second).getRT() - rt_2) < std::abs(rt_earlier - rt_2)))
            {
              rt_earlier = (satellite_it_2->second).getRT();
              intensity_earlier = (satellite_it_2->second).getIntensity();
            }
            
            if ((satellite_it_2->second).getRT() >= rt_2 && (std::abs((satellite_it_2->second).getRT() - rt_2) < std::abs(rt_later - rt_2)))
            {
              rt_later = (satellite_it_2->second).getRT();
              intensity_later = (satellite_it_2->second).getIntensity();
            }
            
          }
          
          // Our target lies on or between two satellites of the 'other' peptide.
          if ((rt_earlier > 0) && (rt_later > 0))
          {
            // linearly interpolated intensity
            double intensity_other;
            if ((rt_2 == rt_earlier) || (rt_later == rt_earlier))
            {
              intensity_other = intensity_earlier;
            }
            else
            {
              intensity_other = intensity_earlier + (intensity_later - intensity_earlier)*(rt_2 - rt_earlier)/(rt_later - rt_earlier);
            }
            
            intensities_light.push_back((satellite_it_1->second).getIntensity());
            intensities_other.push_back(intensity_other);
          }
        }
        
      }
      
      // If less than three matches are found, we cannot reliably calculate the intensity ratio (aka slope) and report the uncorrected intensities.
      if ((intensities_light.size() < 3) || (intensities_other.size() < 3))
      {
        return intensity_peptide;
      }
      
      // determine ratios through linear regression of all corresponding intensities
      Math::LinearRegressionWithoutIntercept linreg;
      linreg.addData(intensities_light, intensities_other);
      double slope = linreg.getSlope();
      
      ratio_peptide.push_back(slope);
    }
    
    // correct peptide intensities
    // The peptide ratios are calculated as linear regression of (spline-interpolated) profile intensities, @see linreg
    // The individual peptide intensities are the sum of the same profile intensities. But the quotient of these peptide intensities
    // is not necessarily the same as the independently calculated ratio from the linear regression. Since the peptide ratio
    // from linear regression is the more accurate one, we correct the two peptide intensities by projecting them onto the ratio.
    // In the end, both peptide ratio from linear regression and the quotient of the peptide intensities are identical.
    std::vector<double> intensity_peptide_corrected;
    if (intensity_peptide.size() == 2)
    {
      double intensity1 = (intensity_peptide[0] + ratio_peptide[1] * intensity_peptide[1]) / (1 + ratio_peptide[1] * ratio_peptide[1]);
      double intensity2 = ratio_peptide[1] * intensity1;
      intensity_peptide_corrected.push_back(intensity1);
      intensity_peptide_corrected.push_back(intensity2);
    }
    else if (intensity_peptide.size() > 2)
    {
      // Now with n instead of two peptide intensities, one needs to project the peptide intensities onto the hyperplane defined
      // by the set of all peptide ratios (TODO). Instead, it is simpler to keep the lightest peptide intensity fixed, and correct
      // only the remaining ones. The correct peptide ratio (from linear regression) is reported on both cases.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
      for (unsigned i = 1; i < intensity_peptide.size(); ++i)
      {
        intensity_peptide_corrected.push_back(ratio_peptide[i] * intensity_peptide[0]);
      }
    }
    else
    {
      // For simple feature detection (singlets) the intensities remain unchanged.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
    }
    
    return intensity_peptide_corrected;
  }

  void FeatureFinderMultiplexAlgorithm::generateMapsCentroided_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, std::vector<std::map<int, GridBasedCluster> >& cluster_results)
  {
    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      // loop over clusters
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();
        
        // Construct a satellite set for the complete peptide multiplet
        // Make sure there are no duplicates, i.e. the same satellite from different filtered peaks.
        std::multimap<size_t, MultiplexSatelliteCentroided > satellites;
        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          MultiplexFilteredPeak peak = filter_results[pattern].getPeak(*point_it);
          // loop over satellites of the peak
          for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = peak.getSatellites().begin(); satellite_it != peak.getSatellites().end(); ++satellite_it)
          {
            // check if this satellite (i.e. these indices) are already in the set
            bool satellite_in_set = false;
            for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_2 = satellites.begin(); satellite_it_2 != satellites.end(); ++satellite_it_2)
            {
              if ((satellite_it_2->second.getRTidx() == satellite_it->second.getRTidx()) && (satellite_it_2->second.getMZidx() == satellite_it->second.getMZidx()))
              {
                satellite_in_set = true;
                break;
              }
            }
            if (satellite_in_set)
            {
              break;
            }
            
            satellites.insert(std::make_pair(satellite_it->first, MultiplexSatelliteCentroided(satellite_it->second.getRTidx(), satellite_it->second.getMZidx())));
          }
        }
        
        // determine peptide intensities
        std::vector<double> peptide_intensities = determinePeptideIntensitiesCentroided_(patterns[pattern], satellites);
        
        // If no reliable peptide intensity can be determined, we do not report the peptide multiplet.
        if (peptide_intensities[0] == -1)
        {
          continue;
        }
        
        std::vector<Feature> features;
        ConsensusFeature consensus;
        bool abort = false;
        
        // construct the feature and consensus maps
        // loop over peptides
        for (size_t peptide = 0; (peptide < patterns[pattern].getMassShiftCount() && !abort); ++peptide)
        {
          // coordinates of the peptide feature
          // RT is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          // m/z is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          Feature feature;
          double rt(0);
          double mz(0);
          double intensity_sum(0);
          
          // loop over isotopes i.e. mass traces of the peptide
          for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
          {
            // find satellites for this isotope i.e. mass trace
            size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
            std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope;
            satellites_isotope = satellites.equal_range(idx);
            
            DBoundingBox<2> mass_trace;
            
            // loop over satellites for this isotope i.e. mass trace
            for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
            {
              // find indices of the peak
              size_t rt_idx = (satellite_it->second).getRTidx();
              size_t mz_idx = (satellite_it->second).getMZidx();
              
              // find peak itself
              MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
              std::advance(it_rt, rt_idx);
              MSSpectrum::ConstIterator it_mz = it_rt->begin();
              std::advance(it_mz, mz_idx);
              
              if (isotope == 0)
              {
                rt += it_rt->getRT() * it_mz->getIntensity();
                mz += it_mz->getMZ() * it_mz->getIntensity();
                intensity_sum += it_mz->getIntensity();
              }
              
              mass_trace.enlarge(it_rt->getRT(), it_mz->getMZ());
            }
            
            if ((mass_trace.width() == 0) || (mass_trace.height() == 0))
            {
              // The mass trace contains only a single point. Add a small margin around
              // the point, otherwise the mass trace is considered empty and not drawn.
              // TODO: Remove the magic number for the margin.
              mass_trace.enlarge(mass_trace.minX() - 0.01, mass_trace.minY() - 0.01);
              mass_trace.enlarge(mass_trace.maxX() + 0.01, mass_trace.maxY() + 0.01);
            }
            
            if (!(mass_trace.isEmpty()))
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.maxY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }
          
          rt /= intensity_sum;
          mz /= intensity_sum;
          
          feature.setRT(rt);
          feature.setMZ(mz);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1.0);
          
          // Check that the feature eluted long enough.
          // DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();    // convex hull of the entire peptide feature
          DBoundingBox<2> box = feature.getConvexHulls()[0].getBoundingBox();    // convex hull of the mono-isotopic mass trace
          if (box.maxX() - box.minX() < static_cast<double>(param_.getValue("algorithm:rt_min")))
          {
            abort = true;
            break;
          }
          
          features.push_back(feature);
          
          if (peptide == 0)
          {
            // The first/lightest peptide acts as anchor of the peptide multiplet consensus.
            // All peptide feature handles are connected to this point.
            consensus.setRT(rt);
            consensus.setMZ(mz);
            consensus.setIntensity(peptide_intensities[peptide]);
            consensus.setCharge(patterns[pattern].getCharge());
            consensus.setQuality(1.0);
          }
          
          FeatureHandle feature_handle;
          feature_handle.setRT(rt);
          feature_handle.setMZ(mz);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus.insert(feature_handle);
          consensus_map_.getColumnHeaders()[peptide].size++;
        }
        
        if (!abort)
        {
          consensus_map_.push_back(consensus);
          for (std::vector<Feature>::iterator it = features.begin(); it != features.end(); ++it)
          {
            feature_map_.push_back(*it);
          }
        }
        
      }
      
    }
    
  }

  void FeatureFinderMultiplexAlgorithm::generateMapsProfile_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, const std::vector<std::map<int, GridBasedCluster> >& cluster_results)
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns.size(), "constructing maps");
    
    
    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      setProgress(++progress);
      
      // loop over clusters
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();
        
        // Construct a satellite set for the complete peptide multiplet
        // Make sure there are no duplicates, i.e. the same satellite from different filtered peaks.
        std::multimap<size_t, MultiplexSatelliteProfile > satellites;
        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          MultiplexFilteredPeak peak = filter_results[pattern].getPeak(*point_it);
          // loop over satellites of the peak
          for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = peak.getSatellitesProfile().begin(); satellite_it != peak.getSatellitesProfile().end(); ++satellite_it)
          {
            satellites.insert(std::make_pair(satellite_it->first, MultiplexSatelliteProfile(satellite_it->second.getRT(), satellite_it->second.getMZ(), satellite_it->second.getIntensity())));
          }
        }
        
        // determine peptide intensities
        std::vector<double> peptide_intensities = determinePeptideIntensitiesProfile_(patterns[pattern], satellites);
        
        // If no reliable peptide intensity can be determined, we do not report the peptide multiplet.
        if (peptide_intensities[0] == -1)
        {
          continue;
        }
        
        std::vector<Feature> features;
        ConsensusFeature consensus;
        bool abort = false;
        
        // construct the feature and consensus maps
        // loop over peptides
        for (size_t peptide = 0; (peptide < patterns[pattern].getMassShiftCount() && !abort); ++peptide)
        {
          // coordinates of the peptide feature
          // RT is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          // m/z is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          Feature feature;
          double rt(0);
          double mz(0);
          double intensity_sum(0);
          
          // loop over isotopes i.e. mass traces of the peptide
          for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
          {
            // find satellites for this isotope i.e. mass trace
            size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
            std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope;
            satellites_isotope = satellites.equal_range(idx);
            
            DBoundingBox<2> mass_trace;
            
            // loop over satellites for this isotope i.e. mass trace
            for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
            {
              if (isotope == 0)
              {
                // Satellites of zero intensity makes sense (borders of peaks), but mess up feature/consensus construction.
                double intensity_temp = (satellite_it->second).getIntensity() + 0.0001;
                
                rt += (satellite_it->second).getRT() * intensity_temp;
                mz += (satellite_it->second).getMZ() * intensity_temp;
                intensity_sum += intensity_temp;
              }
              
              mass_trace.enlarge((satellite_it->second).getRT(), (satellite_it->second).getMZ());
            }
            
            if ((mass_trace.width() == 0) || (mass_trace.height() == 0))
            {
              // The mass trace contains only a single point. Add a small margin around
              // the point, otherwise the mass trace is considered empty and not drawn.
              // TODO: Remove the magic number for the margin.
              mass_trace.enlarge(mass_trace.minX() - 0.01, mass_trace.minY() - 0.01);
              mass_trace.enlarge(mass_trace.maxX() + 0.01, mass_trace.maxY() + 0.01);
            }
            
            if (!(mass_trace.isEmpty()))
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.maxY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }
          
          rt /= intensity_sum;
          mz /= intensity_sum;
          
          feature.setRT(rt);
          feature.setMZ(mz);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1.0);
          
          // Check that the feature eluted long enough.
          // DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();    // convex hull of the entire peptide feature
          DBoundingBox<2> box = feature.getConvexHulls()[0].getBoundingBox();    // convex hull of the mono-isotopic mass trace
          if (box.maxX() - box.minX() < static_cast<double>(param_.getValue("algorithm:rt_min")))
          {
            abort = true;
            break;
          }
          
          features.push_back(feature);
          
          if (peptide == 0)
          {
            // The first/lightest peptide acts as anchor of the peptide multiplet consensus.
            // All peptide feature handles are connected to this point.
            consensus.setRT(rt);
            consensus.setMZ(mz);
            consensus.setIntensity(peptide_intensities[peptide]);
            consensus.setCharge(patterns[pattern].getCharge());
            consensus.setQuality(1.0);
          }
          
          FeatureHandle feature_handle;
          feature_handle.setRT(rt);
          feature_handle.setMZ(mz);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus.insert(feature_handle);
          consensus_map_.getColumnHeaders()[peptide].size++;
        }
        
        if (!abort)
        {
          consensus_map_.push_back(consensus);
          for (std::vector<Feature>::iterator it = features.begin(); it != features.end(); ++it)
          {
            feature_map_.push_back(*it);
          }
        }
        
      }
      
    }
    
    endProgress();
  }

  void FeatureFinderMultiplexAlgorithm::run(MSExperiment& exp, bool progress)
  {
    // parameter section: algorithm, get selected charge range
    String charge_string = param_.getValue("algorithm:charge");
    charge_min_ = charge_string.prefix(':').toInt();
    charge_max_ = charge_string.suffix(':').toInt();
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }
    
    // parameter section: algorithm, get isotopes per peptide range
    String isotopes_per_peptide_string = param_.getValue("algorithm:isotopes_per_peptide");
    isotopes_per_peptide_min_ = isotopes_per_peptide_string.prefix(':').toInt();
    isotopes_per_peptide_max_ = isotopes_per_peptide_string.suffix(':').toInt();
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    progress_ = progress;
    
    // check for empty experimental data
    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }
    
    // update m/z and RT ranges
    exp.updateRanges();
    
    // sort according to RT and MZ
    exp.sortSpectra();
    
    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();
    if (spectrum_type == SpectrumSettings::UNKNOWN)
    {
      spectrum_type = PeakTypeEstimator().estimateType(exp[0].begin(), exp[0].end());
    }

    bool centroided;
    if (param_.getValue("algorithm:spectrum_type")=="automatic")
    {
      centroided = spectrum_type == SpectrumSettings::CENTROID;
    }
    else if (param_.getValue("algorithm:spectrum_type")=="centroid")
    {
      centroided = true;
    }
    else  // "profile"
    {
      centroided = false;
    }
    
    // store experiment in member varaibles
    if (centroided)
    {
      exp.swap(exp_centroid_);
      // exp_profile_ will never be used.
    }
    else
    {
      exp.swap(exp_profile_);
      // exp_centroid_ will be constructed later on.
    }

    /**
     * pick peaks (if input data are in profile mode)
     */
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c; // peak boundaries for chromatograms
    
    if (!centroided)
    {
      PeakPickerHiRes picker;
      Param param = picker.getParameters();
      picker.setLogType(getLogType());
      param.setValue("ms_levels", ListUtils::create<Int>("1"));
      param.setValue("signal_to_noise", 0.0); // signal-to-noise estimation switched off
      picker.setParameters(param);
      
      picker.pickExperiment(exp_profile_, exp_centroid_, boundaries_exp_s, boundaries_exp_c);
    }

    /**
     * filter for peak patterns
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(param_.getValue("algorithm:labels"), param_.getValue("algorithm:missed_cleavages"), label_mass_shift_);
    if (param_.getValue("algorithm:knock_out") == "true")
    {
      generator.generateKnockoutDeltaMasses();
    }
    generator.printSamplesLabelsList();
    generator.printDeltaMassesList();

    std::vector<MultiplexDeltaMasses> masses = generator.getDeltaMassesList();
    std::vector<MultiplexIsotopicPeakPattern> patterns = generatePeakPatterns_(charge_min_, charge_max_, isotopes_per_peptide_max_, masses);
    
    std::vector<MultiplexFilteredMSExperiment> filter_results;
    if (centroided)
    {
      // centroided data
      MultiplexFilteringCentroided filtering(exp_centroid_, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, param_.getValue("algorithm:intensity_cutoff"), param_.getValue("algorithm:rt_band"), param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:peptide_similarity"), param_.getValue("algorithm:averagine_similarity"), param_.getValue("algorithm:averagine_similarity_scaling"), param_.getValue("algorithm:averagine_type"));
      filtering.setLogType(getLogType());
      filter_results = filtering.filter();
    }
    else
    {
      // profile data
      MultiplexFilteringProfile filtering(exp_profile_, exp_centroid_, boundaries_exp_s, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, param_.getValue("algorithm:intensity_cutoff"), param_.getValue("algorithm:rt_band"), param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:peptide_similarity"), param_.getValue("algorithm:averagine_similarity"), param_.getValue("algorithm:averagine_similarity_scaling"), param_.getValue("algorithm:averagine_type"));
      filtering.setLogType(getLogType());
      filter_results = filtering.filter();
    }

    /**
     * cluster filter results
     */
    std::vector<std::map<int, GridBasedCluster> > cluster_results;
    if (centroided)
    {
      // centroided data
      MultiplexClustering clustering(exp_centroid_, param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:rt_typical"), static_cast<double>(param_.getValue("algorithm:rt_min")));
      clustering.setLogType(getLogType());
      cluster_results = clustering.cluster(filter_results);
    }
    else
    {
      // profile data
      MultiplexClustering clustering(exp_profile_, exp_centroid_, boundaries_exp_s, param_.getValue("algorithm:rt_typical"), static_cast<double>(param_.getValue("algorithm:rt_min")));
      clustering.setLogType(getLogType());
      cluster_results = clustering.cluster(filter_results);
    }

    /**
     * construct feature and consensus maps i.e. the final results
     */
    if (centroided)
    {
      //consensus_map.setPrimaryMSRunPath(exp_centroid_.getPrimaryMSRunPath());
      //feature_map.setPrimaryMSRunPath(exp_centroid_.getPrimaryMSRunPath());
      generateMapsCentroided_(patterns, filter_results, cluster_results);
    }
    else
    {
      //consensus_map.setPrimaryMSRunPath(exp_profile_.getPrimaryMSRunPath());
      //feature_map.setPrimaryMSRunPath(exp_profile_.getPrimaryMSRunPath());
      generateMapsProfile_(patterns, filter_results, cluster_results);
    }

  }
  
  FeatureMap& FeatureFinderMultiplexAlgorithm::getFeatureMap()
  {
    return(feature_map_);
  }
  
  ConsensusMap& FeatureFinderMultiplexAlgorithm::getConsensusMap()
  {
    return(consensus_map_);
  }
}
