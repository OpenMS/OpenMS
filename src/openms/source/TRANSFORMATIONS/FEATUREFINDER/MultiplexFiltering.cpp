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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
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

  MultiplexFiltering::MultiplexFiltering(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averigine_type) :
    patterns_(patterns), isotopes_per_peptide_min_(isotopes_per_peptide_min), isotopes_per_peptide_max_(isotopes_per_peptide_max), intensity_cutoff_(intensity_cutoff), rt_band_(rt_band), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit), peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity), averagine_similarity_scaling_(averagine_similarity_scaling), averagine_type_(averigine_type)
  {
    // initialise experiment exp_picked_
    // Any peaks below the intensity cutoff cannot be relevant and are therefore removed.
    // loop over spectra
    for (MSExperiment::ConstIterator it_rt = exp_picked.begin(); it_rt < exp_picked.end(); ++it_rt)
    {
      MSSpectrum<Peak1D> spectrum;
      spectrum.setRT(it_rt->getRT());
      // loop over m/z
      for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        if (it_mz->getIntensity() > intensity_cutoff_)
        {
          Peak1D peak;
          peak.setMZ(it_mz->getMZ());
          peak.setIntensity(it_mz->getIntensity());
          spectrum.push_back(peak);
        }
      }
      exp_picked_.addSpectrum(spectrum);
    }
    exp_picked_.updateRanges();
    
    // initialise blacklist blacklist_
    blacklist_.reserve(exp_picked_.getNrSpectra());
    // loop over spectra
    for (MSExperiment::ConstIterator it_rt = exp_picked_.begin(); it_rt < exp_picked_.end(); ++it_rt)
    {
      std::vector<BlacklistEntry> blacklist_spectrum;
      blacklist_spectrum.reserve(it_rt->size());
      // loop over m/z
      for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        BlacklistEntry entry = white;
        blacklist_spectrum.push_back(entry);
      }
      blacklist_.push_back(blacklist_spectrum);
    }    
  }

  MSExperiment MultiplexFiltering::getWhiteMSExperiment_(White2Original& mapping)
  {
    MSExperiment exp_picked_white;
    // loop over spectra
    for (MSExperiment::ConstIterator it_rt = exp_picked_.begin(); it_rt < exp_picked_.end(); ++it_rt)
    {
      MSSpectrum<Peak1D> spectrum_picked_white;
      spectrum_picked_white.setRT(it_rt->getRT());
      
      std::map<int, int> mapping_spectrum;
      int count = 0;
      // loop over m/z
      for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        if (blacklist_[it_rt - exp_picked_.begin()][it_mz - it_rt->begin()] == white)
        {
          Peak1D peak;
          peak.setMZ(it_mz->getMZ());
          peak.setIntensity(it_mz->getIntensity());
          spectrum_picked_white.push_back(peak);
          
          mapping_spectrum[count] = it_mz - it_rt->begin();
          ++count;
        }
      }
      exp_picked_white.addSpectrum(spectrum_picked_white);
      mapping.push_back(mapping_spectrum);
    }
    exp_picked_white.updateRanges();
    
    return exp_picked_white;
  }
  
  bool MultiplexFiltering::checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak) const
  {
    // Any peak with an intensity greater than <threshold>*<intensity_first_peak> is significant.
    double threshold = 0.3;
    
    // Check that there is a peak.
    int mz_idx = it_rt->findNearest(mz, mz_tolerance);
    if (mz_idx != -1)
    {
      MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
      std::advance(it_mz, mz_idx);
      double intensity = it_mz->getIntensity();
      
      // Check that the peak is significant.
      if (intensity > threshold * intensity_first_peak)
      {
        // There is a high-intensity peak at the position mz.
        return true;
      }
    }

    return false;
  }
  
  bool MultiplexFiltering::filterPeakPositions_(const MSSpectrum<Peak1D>::ConstIterator& it_mz, White2Original& index_mapping, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak) const
  {
    // check if peak position is blacklisted
    if (blacklist_[peak.getRTidx()][peak.getMZidx()] == black)
    {
      return false;
    }
    
    // determine absolute m/z tolerance in Th
    double mz_tolerance;
    if (mz_tolerance_unit_)
    {
      // m/z tolerance in ppm
      // Note that the absolute tolerance varies minimally within an m/z pattern.
      // Hence we calculate it only once here.
      mz_tolerance = it_mz->getMZ() * mz_tolerance_ * 1e-6;
    }
    else
    {
      // m/z tolerance in Th
      mz_tolerance = mz_tolerance_;
    }
    
    // debug output variables
    /*int debug_charge = 4;
    int debug_rt_idx = 35;
    int debug_mz_idx = 6;
    bool debug_now = ((pattern.getCharge() == debug_charge) && (peak.getRTidx() == debug_rt_idx) && (peak.getMZidx() == debug_mz_idx));*/
    
    // debug output
    /*if (debug_now)
    {
      std::cout << "    Inside the peak position filter.\n";
    }*/
    
    // The mass traces of the peptide(s) form a m/z shift pattern. Starting with the mono-isotopic mass trace of each peptide,
    // how long is the series of m/z shifts until the first expected mass trace is missing? We want to see
    // at least isotopes_per_peptide_min_ of these m/z shifts in each peptide. Note that we need to demand subsequent(!) mass traces
    // to be present. Otherwise it would be easy to mistake say a 2+ peptide for a 4+ peptide.
    
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      size_t length = 0;
      bool interrupted = false;
      
      // loop over isotopes i.e. mass traces within the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // calculate m/z shift index in pattern
        size_t idx_mz_shift = peptide * isotopes_per_peptide_max_ + isotope;
        
        double mz_shift = pattern.getMZShiftAt(idx_mz_shift);
        bool found = false;
        
        // loop over spectra in RT band
        for (MSExperiment::ConstIterator it_rt = it_rt_band_begin; it_rt < it_rt_band_end; ++it_rt)
        {
          int i = it_rt->findNearest(it_mz->getMZ() + mz_shift, mz_tolerance);
          
          if (i != -1)
          {
            // Note that unlike primary peaks, satellite peaks are not restricted by the blacklist.
            peak.addSatellite(it_rt - it_rt_begin, index_mapping[it_rt - it_rt_begin][i], idx_mz_shift);
            found = true;
          }          
        }
                
        if (found && (!interrupted))
        {
          ++length;
        }
        else
        {
          interrupted = true;
          if (length < isotopes_per_peptide_min_)
          {
            return false;
          }
        }
      }
    }
    
    // Check that there is no significant peak (aka zeroth peak) to the left of the mono-isotopic peak (aka first peak).
    // Further check that there is no mistaken charge state identity. For example, check that a 2+ pattern isn't really a 4+ or 6+ pattern.
    // Let's use the double m/z tolerance when checking for these peaks.
    
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      MSExperiment::ConstIterator it_rt = it_rt_begin;
      std::advance(it_rt, peak.getRTidx());

      // Check that there is a first i.e. mono-isotopic peak for this peptide.
      double mz_first_peak = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_);
      int mz_idx_first_peak = it_rt->findNearest(mz_first_peak, mz_tolerance);
      if (mz_idx_first_peak != -1)
      {
        MSSpectrum<Peak1D>::ConstIterator it_mz_first_peak = it_rt->begin();
        std::advance(it_mz_first_peak, mz_idx_first_peak);
        double intensity_first_peak = it_mz_first_peak->getIntensity();

        double mz;
        
        // Check that there is a zeroth peak.
        mz = peak.getMZ() + 2 * pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_) - pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1);        
        if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak))
        {
          return false;
        }
        
        // Check mistaken charge state identities
        // We are searching the patterns in the order of the most common occurrence (and not decreasing charge state). 
        // That can lead to mistaken charge state identities. Here we check that this is not the case.

        if (pattern.getCharge() == 2)
        {          
          // Is the 2+ pattern really a 4+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)/2 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/2;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak))
          {
            return false;
          }
          
          // Is the 2+ pattern really a 6+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)*2/3 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/3;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak))
          {
            return false;
          }
        }
        
        if (pattern.getCharge() == 3)
        {
          // Is the 3+ pattern really a 6+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)/2 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/2;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak))
          {
            return false;
          }
        }
        
        if (pattern.getCharge() == 1)
        {
          for (int c = 2; c < 7; ++c)
          {
            // Is the 1+ pattern really a c+ pattern?
            mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)*(c-1)/c + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/c;
            if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak))
            {
              return false;
            }
          }
        }
        
      }
            
    }
    
    // Automatically length >= isotopes_per_peptide_min_
    return true;
  }

  void MultiplexFiltering::blacklistPeak_(const MultiplexFilteredPeak& peak)
  {
    // blacklist satellite peaks
    // Note that the primary peak is one of them.
    std::multimap<size_t, MultiplexSatellite > satellites = peak.getSatellites();
    for (std::multimap<size_t, MultiplexSatellite >::const_iterator it = satellites.begin(); it != satellites.end(); ++it)
    {
      if (it->first == 0)    // first = m/z shift index in <MultiplexIsotopicPeakPattern>
      {
        // Any mono-isotopic peaks of the lightest peptide can pass the filter in this pattern.
        blacklist_[(it->second).getRTidx()][(it->second).getMZidx()] = grey;    // first = RT index in original <MSExperiment>, second = m/z index in spectrum
      }
      else
      {
        blacklist_[(it->second).getRTidx()][(it->second).getMZidx()] = black;
      }
    }
  }
  
  void MultiplexFiltering::blacklistPeak2_(const MultiplexFilteredPeak& peak, unsigned pattern_idx)
  {
    // determine absolute m/z tolerance in Th
    double mz_tolerance;
    if (mz_tolerance_unit_)
    {
      // m/z tolerance in ppm
      // Note that the absolute tolerance varies minimally within an m/z pattern.
      // Hence we calculate it only once here.
      mz_tolerance = peak.getMZ() * mz_tolerance_ * 1e-6;
    }
    else
    {
      // m/z tolerance in Th
      mz_tolerance = mz_tolerance_;
    }    

    // Determine the boundaries for each of the mass traces.
    std::multimap<size_t, MultiplexSatellite > satellites = peak.getSatellites();    
    // <rt_boundaries> is a map from the mass trace index to the spectrum indices for beginning and end of the mass trace.
    std::map<size_t, std::pair<size_t, size_t> > rt_boundaries;
    // loop over satellites
    for (std::multimap<size_t, MultiplexSatellite >::const_iterator it = satellites.begin(); it != satellites.end(); ++it)
    {
      size_t idx_masstrace = it->first;    // mass trace index i.e. the index within the peptide multiplet pattern
      if (rt_boundaries.find(idx_masstrace) == rt_boundaries.end())
      {
        // That's the first satellite within this mass trace.
        rt_boundaries[idx_masstrace] = std::make_pair((it->second).getRTidx(), (it->second).getRTidx());
      }
      else
      {
        // We have seen a satellite of this mass trace before.
        size_t idx_min = std::min((it->second).getRTidx(), rt_boundaries[idx_masstrace].first);
        size_t idx_max = std::max((it->second).getRTidx(), rt_boundaries[idx_masstrace].second);
        
        rt_boundaries[idx_masstrace] = std::make_pair(idx_min, idx_max);
      }      
    }
        
    // Blacklist all peaks along the mass traces
    // loop over mass traces (i.e. the mass trace boundaries)
    for (std::map<size_t, std::pair<size_t, size_t> >::const_iterator it = rt_boundaries.begin(); it != rt_boundaries.end(); ++it)
    {
      double mz = peak.getMZ() + patterns_[pattern_idx].getMZShiftAt(it->first);
      
      MSExperiment::ConstIterator it_rt_begin = exp_picked_.begin() + (it->second).first;
      it_rt_begin = exp_picked_.RTBegin(it_rt_begin->getRT() - rt_band_);
      
      MSExperiment::ConstIterator it_rt_end = exp_picked_.begin() + (it->second).second;
      it_rt_end = exp_picked_.RTBegin(it_rt_end->getRT() + rt_band_);
      
      // prepare for loop
      if (it_rt_end != exp_picked_.end())
      {
        ++it_rt_end;
      }

      // loop over RT along the mass trace
      for (MSExperiment::ConstIterator it_rt = it_rt_begin; it_rt < it_rt_end; ++it_rt)
      {        
        int idx_mz = it_rt->findNearest(mz, mz_tolerance);
       
        if (idx_mz != -1)
        {
          if (it->first == 0)
          {
            // Any mono-isotopic peaks of the lightest peptide can pass the filter in this pattern.
            blacklist_[it_rt - exp_picked_.begin()][idx_mz] = grey;
          }
          else
          {
            blacklist_[it_rt - exp_picked_.begin()][idx_mz] = black;
          }
        }
      }
    
    }
    
  }
  
  void MultiplexFiltering::ungreyBlacklist_()
  {
    // loop over spectra
    for (MSExperiment::ConstIterator it_rt = exp_picked_.begin(); it_rt < exp_picked_.end(); ++it_rt)
    {
      // loop over m/z
      for (MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
      {
        if (blacklist_[it_rt - exp_picked_.begin()][it_mz - it_rt->begin()] == grey)
        {
          blacklist_[it_rt - exp_picked_.begin()][it_mz - it_rt->begin()] = black;
        }
      }
    }
  }
  
  bool MultiplexFiltering::filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const
  {
    // construct averagine distribution
    // Note that the peptide(s) are very close in mass. We therefore calculate the averagine distribution only once (for the lightest peptide).
    double mass = peak.getMZ() * pattern.getCharge();
    IsotopeDistribution distribution;
    distribution.setMaxIsotope(isotopes_per_peptide_max_);
    if (averagine_type_ == "peptide")
    {
        distribution.estimateFromPeptideWeight(mass);
    }
    else if (averagine_type_ == "RNA")
    {
        distribution.estimateFromRNAWeight(mass);
    }
    else if (averagine_type_ == "DNA")
    {
        distribution.estimateFromDNAWeight(mass);
    }
    else
    {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid averagine type.");
    }
    
    // debug output variables
    /*int debug_charge = 2;
    size_t debug_rt_idx = 39;
    size_t debug_mz_idx = 130;
    bool debug_now = ((pattern.getCharge() == debug_charge) && (peak.getRTidx() == debug_rt_idx) && (peak.getMZidx() == debug_mz_idx));*/
    
    // debug output
    /*if (debug_now)
    {
      std::cout << "Inside the Averagine Filter.\n";
    }*/
    
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
        std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites;
        satellites = peak.getSatellites().equal_range(idx);
              
        int count = 0;
        double sum_intensities = 0;
        
        // loop over satellites in mass trace
        for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it = satellites.first; satellite_it != satellites.second; ++satellite_it)
        {
          // find indices of the peak
          size_t rt_idx = (satellite_it->second).getRTidx();
          size_t mz_idx = (satellite_it->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt = exp_picked_.begin();
          std::advance(it_rt, rt_idx);
          MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
          std::advance(it_mz, mz_idx);
                    
          ++count;
          sum_intensities += it_mz->getIntensity();
        }
        
        if (count > 0)
        {
          intensities_model.push_back(distribution.getContainer()[isotope].second);
          intensities_data.push_back(sum_intensities/count);
        }
        
        // debug output
        /*if (debug_now)
        {
          std::cout << "    peptide = " << peptide << "    isotope = " << isotope << "    count = " << count << "    average intensity = " << sum_intensities/count << "    averagine intensity = " << distribution.getContainer()[isotope].second << "\n";
        }*/
                
      }
      
      // Calculate Pearson and Spearman rank correlations
      if ((intensities_model.size() < isotopes_per_peptide_min_) || (intensities_data.size() < isotopes_per_peptide_min_))
      {
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
      }
      double correlation_Pearson = OpenMS::Math::pearsonCorrelationCoefficient(intensities_model.begin(), intensities_model.end(), intensities_data.begin(), intensities_data.end());
      double correlation_Spearman = OpenMS::Math::rankCorrelationCoefficient(intensities_model.begin(), intensities_model.end(), intensities_data.begin(), intensities_data.end());

      // debug output
      /*if (debug_now)
      {
        std::cout << "        Pearson correlation = " << correlation_Pearson << "    rank correlation = " << correlation_Spearman << "\n";
      }*/
      
      if ((correlation_Pearson < averagine_similarity_) || (correlation_Spearman < averagine_similarity_))
      {
        return false;
      }
    }
    
    return true;
  }

  bool MultiplexFiltering::filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const
  {
    if (pattern.getMassShiftCount() < 2)
    {
      // filter irrelevant for singlet feature detection
      return true;
    }
    
    // debug output variables
    int debug_charge = 4;
    size_t debug_rt_idx = 35;
    size_t debug_mz_idx = 6;
    bool debug_now = ((pattern.getCharge() == debug_charge) && (peak.getRTidx() == debug_rt_idx) && (peak.getMZidx() == debug_mz_idx));
    
    // debug output
    /*if (debug_now)
    {
      std::cout << "Inside the Peptide Correlation Filter.\n";
    }*/
    
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
                    
          std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_1;
          std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_2;
          satellites_1 = peak.getSatellites().equal_range(idx_1);
          satellites_2 = peak.getSatellites().equal_range(idx_2);
          
          //std::cout << "isotope = " << isotope << "    idx_1 = " << idx_1 << "    idx_2 = " << idx_2 << "\n";
          
          // loop over satellites in mass trace 1
          for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_1 = satellites_1.first; satellite_it_1 != satellites_1.second; ++satellite_it_1)
          {
            size_t rt_idx_1 = (satellite_it_1->second).getRTidx();
  
            //std::cout << "    Inside Loop.    rt_idx_1 = " << rt_idx_1 << "\n";
            
            // loop over satellites in mass trace 2
            for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_2 = satellites_2.first; satellite_it_2 != satellites_2.second; ++satellite_it_2)
            {
              size_t rt_idx_2 = (satellite_it_2->second).getRTidx();
  
              //std::cout << "        Inside Second Loop.    rt_idx_2 = " << rt_idx_2 << "\n";

              if (rt_idx_1 == rt_idx_2)
              {
                std::cout << "        Both rt_idx equal.\n";
                
                size_t mz_idx_1 = (satellite_it_1->second).getMZidx();
                size_t mz_idx_2 = (satellite_it_2->second).getMZidx();
                
                // find peak itself
                MSExperiment::ConstIterator it_rt_1 = exp_picked_.begin();
                MSExperiment::ConstIterator it_rt_2 = exp_picked_.begin();
                std::advance(it_rt_1, rt_idx_1);
                std::advance(it_rt_2, rt_idx_2);
                MSSpectrum<Peak1D>::ConstIterator it_mz_1 = it_rt_1->begin();
                MSSpectrum<Peak1D>::ConstIterator it_mz_2 = it_rt_2->begin();
                std::advance(it_mz_1, mz_idx_1);
                std::advance(it_mz_2, mz_idx_2);
  
                intensities_1.push_back(it_mz_1->getIntensity());
                intensities_2.push_back(it_mz_2->getIntensity());
              }
            }
          }
        }
        
        // It is well possible that no corresponding satellite peaks exist, in which case the filter fails.
        if ((intensities_1.size() == 0) || (intensities_2.size() == 0))
        {
          return false;
        }
        
        // calculate correlation between peak insities in peptides 1 and 2
        double correlation_Pearson = OpenMS::Math::pearsonCorrelationCoefficient(intensities_1.begin(), intensities_1.end(), intensities_2.begin(), intensities_2.end());
        double correlation_Spearman = OpenMS::Math::rankCorrelationCoefficient(intensities_1.begin(), intensities_1.end(), intensities_2.begin(), intensities_2.end());

        // debug output
        /*if (debug_now)
        {
          std::cout << "        Pearson correlation = " << correlation_Pearson << "    rank correlation = " << correlation_Spearman << "\n";
          //std::cout << "        Pearson correlation = " << correlation_Pearson << "\n";
        }*/
        
        if ((correlation_Pearson < peptide_similarity_) || (correlation_Spearman < peptide_similarity_))
        //if (correlation_Pearson < peptide_similarity_)
        {
          return false;
        }
      }
    }
    
    return true;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  bool MultiplexFiltering::monoIsotopicPeakIntensityFilter_(const MultiplexIsotopicPeakPattern& pattern, int spectrum_index, const vector<int>& mz_shifts_actual_indices) const
  {
    MSExperiment::ConstIterator it_rt = exp_picked_.begin() + spectrum_index;
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      int peak_index = mz_shifts_actual_indices[peptide * isotopes_per_peptide_max_];
      if (peak_index == -1)
      {
        // peak not found
        return true;
      }
      MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin() + peak_index;
      if (it_mz->getIntensity() < intensity_cutoff_)
      {
        // below intensity threshold
        return true;
      }
    }
    return false;
  }

  bool MultiplexFiltering::zerothPeakFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual) const
  {
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // scaling factor for the zeroth peak intensity
      // (The zeroth peak is problematic if its intensity exceeds zero_scaling * intensity of mono-isotopic peak.)
      double zero_scaling = 0.7;
      if (boost::math::isnan(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1)]))
      {
        // zeroth peak not found
        continue;
      }
      else if (boost::math::isnan(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + 1]))
      {
        // first peak not found
        return true;
      }
      else if (intensities_actual[peptide * (isotopes_per_peptide_max_ + 1)] > zero_scaling * intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + 1])
      {
        return true;
      }
    }

    return false;
  }

  bool MultiplexFiltering::peptideSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline) const
  {
    std::vector<double> isotope_pattern_1;
    std::vector<double> isotope_pattern_2;
    for (unsigned peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
      {
        if (boost::math::isnan(intensities_actual[isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern_1.push_back(0);
        }
        else
        {
          isotope_pattern_1.push_back(intensities_actual[isotope + 1]);
        }
        if (boost::math::isnan(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern_2.push_back(0);
        }
        else
        {
          isotope_pattern_2.push_back(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1]);
        }
      }
      if (getPatternSimilarity_(isotope_pattern_1, isotope_pattern_2) < peptide_similarity_)
      {
        return false;
      }
    }

    return true;
  }

  bool MultiplexFiltering::averagineSimilarityFilter_(const MultiplexIsotopicPeakPattern& pattern, const vector<double>& intensities_actual, int peaks_found_in_all_peptides_spline, double mz) const
  {
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
    
    for (unsigned peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      vector<double> isotope_pattern;
      for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
      {
        if (boost::math::isnan(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1]))
        {
          // no peak found, hence assume the intensity at this position to be zero
          isotope_pattern.push_back(0);
        }
        else
        {
          isotope_pattern.push_back(intensities_actual[peptide * (isotopes_per_peptide_max_ + 1) + isotope + 1]);
        }
      }
      if (getAveragineSimilarity_(isotope_pattern, mz * pattern.getCharge()) < similarity)
      {
        return false;
      }
    }

    return true;
  }

  int MultiplexFiltering::getPeakIndex_(const std::vector<double>& peak_position, int start, double mz, double scaling) const
  {
    const double tolerance_th = mz_tolerance_unit_ ? (scaling * mz_tolerance_ / 1000000) * peak_position[start] : scaling * mz_tolerance_;
    const double mz_min = mz - tolerance_th;
    const double mz_max = mz + tolerance_th;

    std::vector<double>::const_iterator lb = std::lower_bound(peak_position.begin(), peak_position.end(), mz_min);
    std::vector<double>::const_iterator ub = std::upper_bound(lb, peak_position.end(), mz_max);

    double smallest_error = scaling * mz_tolerance_; // initialize to the maximum  allowed error 
    int smallest_error_index = -1;

    for (; lb != ub; ++lb)
    {
      const double error = abs(*lb - mz);
      if (error <= smallest_error)
      {
        smallest_error = error;
        smallest_error_index = lb - peak_position.begin();
      }
    }    

    return smallest_error_index;
  }

  double MultiplexFiltering::getPatternSimilarity_(const vector<double>& pattern1, const vector<double>& pattern2) const
  {
    if (pattern1.empty() || pattern2.empty())
    {
      return std::numeric_limits<double>::quiet_NaN();
    }

    return OpenMS::Math::pearsonCorrelationCoefficient(pattern1.begin(), pattern1.end(), pattern2.begin(), pattern2.end());
  }


  double MultiplexFiltering::getAveragineSimilarity_(const vector<double>& pattern, double m) const

  {
    // construct averagine distribution
    IsotopeDistribution distribution;
    vector<double> averagine_pattern;
    distribution.setMaxIsotope(pattern.size());
    if (averagine_type_ == "peptide")
    {
        distribution.estimateFromPeptideWeight(m);
    }
    else if (averagine_type_ == "RNA")
    {
        distribution.estimateFromRNAWeight(m);
    }
    else if (averagine_type_ == "DNA")
    {
        distribution.estimateFromDNAWeight(m);
    }
    else
    {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Averagine type unrecognized.");;
    }

    for (IsotopeDistribution::Iterator it = distribution.begin(); it != distribution.end(); ++it)
    {
      averagine_pattern.push_back(it->second);
    }

    return getPatternSimilarity_(pattern, averagine_pattern);
  }

}
