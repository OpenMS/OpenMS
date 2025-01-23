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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <sstream>
#include <utility>

//#define DEBUG_FFMULTIPLEX

using namespace std;

namespace OpenMS
{

  MultiplexFilteringProfile::MultiplexFilteringProfile(MSExperiment& exp_profile, const MSExperiment& exp_centroided, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type) :
    MultiplexFiltering(exp_centroided, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, std::move(averagine_type))
  {
    // initialise peak boundaries
    // In the MultiplexFiltering() constructor we initialise the centroided experiment exp_centroided_.
    // (We run a simple intensity filter. Peaks below the intensity cutoff can be discarded right from the start.)
    // Now we still need to discard boundaries of low intensity peaks, in order to preserve the one-to-one mapping between peaks and boundaries.
    boundaries_.reserve(boundaries.size());
    // loop over spectra and boundaries
    for (const auto &it_rt : exp_centroided)
    {
      size_t idx_rt = &it_rt - &exp_centroided[0];
      
      // new boundaries of a single spectrum
      std::vector<PeakPickerHiRes::PeakBoundary> boundaries_temp;
      
      // loop over m/z peaks and boundaries
      for (const auto &it_mz : it_rt)
      {
        size_t idx_mz = &it_mz - &it_rt[0];
        
        if (it_mz.getIntensity() > intensity_cutoff_)
        {
          boundaries_temp.push_back(boundaries[idx_rt][idx_mz]);

          // Check consistency of peaks and their peak boundaries, i.e. check that the peak lies in the boundary interval.
          if (boundaries[idx_rt][idx_mz].mz_min > it_mz.getMZ() || it_mz.getMZ() > boundaries[idx_rt][idx_mz].mz_max)
          {
            throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
          }
        }
      }
      
      boundaries_.push_back(boundaries_temp);
    }
    
    if (exp_profile.size() != exp_centroided.size())
    {
      stringstream stream;
      stream << "Profile and centroided data do not contain same number of spectra. (";
      stream << exp_profile.size();
      stream << "!=";
      stream << exp_centroided.size();
      stream << ")";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
    }

    if (exp_centroided.size() != boundaries.size())
    {
      stringstream stream;
      stream << "Centroided data and the corresponding list of peak boundaries do not contain same number of spectra. (";
      stream << exp_centroided.size();
      stream << "!=";
      stream << boundaries.size();
      stream << ")";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, stream.str());
    }
    
    // spline interpolate the profile data
    for (MSExperiment::Iterator it = exp_profile.begin(); it < exp_profile.end(); ++it)
    {
      exp_spline_profile_.emplace_back(*it);
    }
    
    // TODO: Constructing the navigators here instead in the beginning of the filter() method results in segmentation faults. Why?

  }
  
  vector<MultiplexFilteredMSExperiment> MultiplexFilteringProfile::filter()
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns_.size() * exp_spline_profile_.size(), "filtering LC-MS data");

    // list of filter results for each peak pattern
    std::vector<MultiplexFilteredMSExperiment> filter_results;
    
#ifdef DEBUG_FFMULTIPLEX
    // clock for monitoring run time performance
    unsigned int start = clock();
#endif
    
    // construct navigators for all spline spectra
    std::vector<SplineInterpolatedPeaks::Navigator> navigators;
    for (SplineInterpolatedPeaks& spl : exp_spline_profile_)
    {
      SplineInterpolatedPeaks::Navigator nav = spl.getNavigator();
      navigators.push_back(nav);
    }
    
    // loop over all patterns
    for (unsigned pattern_idx = 0; pattern_idx < patterns_.size(); ++pattern_idx)
    {
      // current pattern
      MultiplexIsotopicPeakPattern pattern = patterns_[pattern_idx];
      
      // data structure storing peaks which pass all filters
      MultiplexFilteredMSExperiment result;

      // update white experiment
      updateWhiteMSExperiment_();
      
      // loop over spectra
      // loop simultaneously over RT in the spline interpolated profile and (white) centroided experiment (including peak boundaries)
      for (const auto &it_rt : exp_centroided_white_)
      {
        // retention time
        double rt = it_rt.getRT();
        // spectral index in exp_centroided_white_, boundaries_ and exp_spline_profile_
        size_t idx_rt = &it_rt - &exp_centroided_white_[0];
        
        // skip empty spectra
        if (it_rt.empty() || boundaries_[idx_rt].empty() || exp_spline_profile_[idx_rt].size() == 0)
        {
          continue;
        }
        
        setProgress(++progress);
        
        MSExperiment::ConstIterator it_rt_picked_band_begin = exp_centroided_white_.RTBegin(rt - rt_band_/2);
        MSExperiment::ConstIterator it_rt_picked_band_end = exp_centroided_white_.RTEnd(rt + rt_band_/2);
        
        // loop over mz
        #pragma omp parallel for
        for (SignedSize s = 0; s < (SignedSize) it_rt.size(); s++)
        {
          double mz = it_rt[s].getMZ();
          MultiplexFilteredPeak peak(mz, rt, exp_centroided_mapping_[idx_rt][s], idx_rt);
          
          if (!(filterPeakPositions_(mz, exp_centroided_white_.begin(), it_rt_picked_band_begin, it_rt_picked_band_end, pattern, peak)))
          {
            continue;
          }
          
          size_t mz_idx = exp_centroided_mapping_[idx_rt][s];
          double peak_min = boundaries_[idx_rt][mz_idx].mz_min;
          double peak_max = boundaries_[idx_rt][mz_idx].mz_max;
          
          //double rt_peak = peak.getRT();
          double mz_peak = peak.getMZ();

          std::multimap<size_t, MultiplexSatelliteCentroided > satellites = peak.getSatellites();
          
          // Arrangement of peaks looks promising. Now scan through the spline fitted profile data around the peak i.e. from peak boundary to peak boundary.
          for (double mz_profile = peak_min; mz_profile < peak_max; mz_profile = navigators[idx_rt].getNextPos(mz_profile))
          {
            // determine m/z shift relative to the centroided peak at which the profile data will be sampled
            double mz_shift = mz_profile - mz_peak;

            std::multimap<size_t, MultiplexSatelliteProfile > satellites_profile;

            // construct the set of spline-interpolated satellites for this specific mz_profile
            for (const auto &satellite_it : satellites)
            {
              // find indices of the peak
              size_t rt_idx = (satellite_it.second).getRTidx();
              size_t mz_idx = (satellite_it.second).getMZidx();
              
              // find peak itself
              MSExperiment::ConstIterator it_rt = exp_centroided_.begin();
              std::advance(it_rt, rt_idx);
              MSSpectrum::ConstIterator it_mz = it_rt->begin();
              std::advance(it_mz, mz_idx);
              
              double rt_satellite = it_rt->getRT();
              double mz_satellite = it_mz->getMZ();
              
              // determine m/z and corresponding intensity
              double mz = mz_satellite + mz_shift;
              double intensity = navigators[rt_idx].eval(mz);
              
              satellites_profile.insert(std::make_pair(satellite_it.first, MultiplexSatelliteProfile(rt_satellite, mz, intensity)));
            }
            
            if (!(filterAveragineModel_(pattern, peak, satellites_profile)))
            {
              continue;
            }
            
            if (!(filterPeptideCorrelation_(pattern, satellites_profile)))
            {
              continue;
            }
            
            /**
             * All filters passed.
             */
            
            // add the satellite data points to the peak
            for (const auto &it : satellites_profile)
            {
              peak.addSatelliteProfile(it.second, it.first);
            }
            
          }
          
          // If some satellite data points passed all filters, we can add the peak to the filter result.
          if (peak.sizeProfile() > 0)
          {
            #pragma omp critical
            {
              result.addPeak(peak);
              blacklistPeak_(peak, pattern_idx);
            };
          }
          
        }
        
      }
 
#ifdef DEBUG_FFMULTIPLEX
      // write filtered peaks to debug output
      std::stringstream debug_out;
      debug_out << "filter_result_" << pattern_idx << ".consensusXML";
      result.writeDebugOutput(exp_centroided_, debug_out.str());
#endif
      
      // add results of this pattern to list
      filter_results.push_back(result);
    }
    
#ifdef DEBUG_FFMULTIPLEX
    // clock for monitoring run time performance
    OPENMS_LOG_INFO << "\nThe filtering step of the algorithm took " << (float)(clock()-start)/CLOCKS_PER_SEC << " seconds.\n\n";
#endif

    endProgress();

    return filter_results;
  }
  
  std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& MultiplexFilteringProfile::getPeakBoundaries()
  {
    return boundaries_;
  }

  bool MultiplexFilteringProfile::filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const
  {
    // construct averagine distribution
    // Note that the peptide(s) are very close in mass. We therefore calculate the averagine distribution only once (for the lightest peptide).
    double mass = peak.getMZ() * pattern.getCharge();
    CoarseIsotopePatternGenerator solver(isotopes_per_peptide_max_);
    IsotopeDistribution distribution;
    if (averagine_type_ == "peptide")
    {
      distribution = solver.estimateFromPeptideWeight(mass);
    }
    else if (averagine_type_ == "RNA")
    {
      distribution = solver.estimateFromRNAWeight(mass);
    }
    else if (averagine_type_ == "DNA")
    {
      distribution = solver.estimateFromDNAWeight(mass);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid averagine type.");
    }   
    
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // intensities for the Pearson and Spearman rank correlations
      std::vector<double> intensities_model;
      std::vector<double> intensities_data;
      
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites;
        satellites = satellites_profile.equal_range(idx);
        
        int count = 0;
        double sum_intensities = 0;
        
        // loop over satellites in mass trace
        for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator it = satellites.first; it != satellites.second; ++it)
        {
          ++count;
          sum_intensities += (it->second).getIntensity();
        }
        
        if (count > 0)
        {
          intensities_model.push_back(distribution[isotope].getIntensity());
          intensities_data.push_back(sum_intensities/count);
        }
        
      }
      
      // Use a more restrictive averagine similarity when we are searching for peptide singlets.
      double similarity;
      if (pattern.getMassShiftCount() == 1)
      {
        // We are detecting peptide singlets.
        similarity = averagine_similarity_ + averagine_similarity_scaling_*(1 - averagine_similarity_);
      }
      else
      {
        // We are detecting peptide doublets or triplets or ...
        similarity = averagine_similarity_;
      }
      
      // Calculate Pearson and Spearman rank correlations
      if ((intensities_model.size() < isotopes_per_peptide_min_) || (intensities_data.size() < isotopes_per_peptide_min_))
      {
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
      }
      double correlation_Pearson = OpenMS::Math::pearsonCorrelationCoefficient(intensities_model.begin(), intensities_model.end(), intensities_data.begin(), intensities_data.end());
      double correlation_Spearman = OpenMS::Math::rankCorrelationCoefficient(intensities_model.begin(), intensities_model.end(), intensities_data.begin(), intensities_data.end());

      if ((correlation_Pearson < similarity) || (correlation_Spearman < similarity))
      {
        return false;
      }

      
    }
    
    return true;
  }
  
  bool MultiplexFilteringProfile::filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern,
                                                            const std::multimap<size_t, MultiplexSatelliteProfile >& satellites_profile) const
  {
    if (pattern.getMassShiftCount() < 2)
    {
      // filter irrelevant for singlet feature detection
      return true;
    }

    // We will calculate the correlations between all possible peptide combinations.
    // For example (light, medium), (light, heavy) and (medium, heavy) in the case of triplets.
    // If one of the correlations is below the <peptide_similarity_> limit, the filter fails.
    
    // loop over the first peptide
    for (size_t peptide_1 = 0; peptide_1 < pattern.getMassShiftCount() - 1; ++peptide_1)
    {
      // loop over the second peptide
      for (size_t peptide_2 = peptide_1 + 1; peptide_2 < pattern.getMassShiftCount(); ++peptide_2)
      {
        std::vector<double> intensities_1;
        std::vector<double> intensities_2;
        
        // loop over isotopes i.e. mass traces of both peptides
        for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
        {
          size_t idx_1 = peptide_1 * isotopes_per_peptide_max_ + isotope;
          size_t idx_2 = peptide_2 * isotopes_per_peptide_max_ + isotope;
          
          auto satellites_1 = satellites_profile.equal_range(idx_1);
          auto satellites_2 = satellites_profile.equal_range(idx_2);
          
          // loop over satellites in mass trace 1
          for (auto satellite_it_1 = satellites_1.first; satellite_it_1 != satellites_1.second; ++satellite_it_1) //OMS_CODING_TEST_EXCLUDE
          {
            double rt_1 = (satellite_it_1->second).getRT();
            
            // loop over satellites in mass trace 2
            for (auto satellite_it_2 = satellites_2.first; satellite_it_2 != satellites_2.second; ++satellite_it_2) //OMS_CODING_TEST_EXCLUDE
            {
              double rt_2 = (satellite_it_2->second).getRT();
              
              if (rt_1 == rt_2)
              {
                intensities_1.push_back((satellite_it_1->second).getIntensity());
                intensities_2.push_back((satellite_it_2->second).getIntensity());
              }

            }
            
          }

        }

        // It is well possible that no corresponding satellite peaks exist, in which case the filter fails.
        if ((intensities_1.empty()) || (intensities_2.empty()))
        {
          return false;
        }
        
        // calculate correlation between peak intensities in peptides 1 and 2
        double correlation_Pearson = OpenMS::Math::pearsonCorrelationCoefficient(intensities_1.begin(), intensities_1.end(), intensities_2.begin(), intensities_2.end());
        double correlation_Spearman = OpenMS::Math::rankCorrelationCoefficient(intensities_1.begin(), intensities_1.end(), intensities_2.begin(), intensities_2.end());
        
        if ((correlation_Pearson < peptide_similarity_) || (correlation_Spearman < peptide_similarity_))
        //if (correlation_Pearson < peptide_similarity_)
        {
          return false;
        }

      }
      
    }
    
    return true;
  }
  
}
