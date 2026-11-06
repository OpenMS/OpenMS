// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  namespace
  {
    /**
      The same as spectrum.findNearest(mz, tolerance), but the search for the first peak at or right of @p mz starts
      from @p cursor, where the previous search on the spectrum ended: peaks are searched at increasing m/z, so it is
      usually close. (Any cursor gives the same result: unless all peaks left of it lie left of @p mz, the whole
      spectrum is searched.) Afterwards, @p cursor is that first peak at or right of @p mz.
    */
    Int findNearestFrom(const MSSpectrum& spectrum, double mz, double tolerance, Size& cursor)
    {
      if (spectrum.empty())
      {
        return -1;
      }
      const Size n = spectrum.size();
      Size lo = std::min(cursor, n);
      MSSpectrum::ConstIterator it;
      if (lo > 0 && !(spectrum[lo - 1].getMZ() < mz))
      {
        it = spectrum.MZBegin(mz);
      }
      else
      {
        // All peaks before lo lie left of mz. Double the step to the right until a peak at or right of mz is passed.
        Size hi = lo;
        for (Size step = 1; hi < n && spectrum[hi].getMZ() < mz; step *= 2)
        {
          lo = hi + 1;
          hi += step;
        }
        it = spectrum.MZBegin(spectrum.begin() + lo, mz, spectrum.begin() + std::min(hi, n));
      }
      cursor = static_cast<Size>(it - spectrum.begin());

      // as MSSpectrum::findNearest(mz) and findNearest(mz, tolerance)
      Size i;
      if (it == spectrum.begin())
      {
        i = 0;
      }
      else if (it == spectrum.end())
      {
        i = n - 1;
      }
      else
      {
        MSSpectrum::ConstIterator it2 = it;
        --it2;
        i = (std::fabs(it->getMZ() - mz) < std::fabs(it2->getMZ() - mz)) ? Size(it - spectrum.begin()) : Size(it2 - spectrum.begin());
      }
      const double found_mz = spectrum[i].getMZ();
      if (found_mz >= mz - tolerance && found_mz <= mz + tolerance)
      {
        return static_cast<Int>(i);
      }
      return -1;
    }
  }

  MultiplexFiltering::MultiplexFiltering(const MSExperiment& exp_centroided, const std::vector<MultiplexIsotopicPeakPattern>& patterns,
                                         int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band,
                                         double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity,
                                         double averagine_similarity_scaling, std::string averagine_type) :
  patterns_(patterns), isotopes_per_peptide_min_(isotopes_per_peptide_min), isotopes_per_peptide_max_(isotopes_per_peptide_max),
  intensity_cutoff_(intensity_cutoff), rt_band_(rt_band), mz_tolerance_(mz_tolerance), mz_tolerance_unit_in_ppm_(mz_tolerance_unit),
  peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity),
  averagine_similarity_scaling_(averagine_similarity_scaling), averagine_type_(std::move(averagine_type))
  {
    // initialise experiment exp_centroided_
    // Any peaks below the intensity cutoff cannot be relevant. They are therefore removed resulting in reduced memory footprint and runtime.
    exp_centroided_.reserve(exp_centroided.getNrSpectra());
    // loop over spectra
    for (const auto &it_rt : exp_centroided)
    {
      MSSpectrum spectrum;
      spectrum.setRT(it_rt.getRT());
      // loop over m/z
      for (const auto &it_mz : it_rt)
      {
        if (it_mz.getIntensity() > intensity_cutoff_)
        {
          spectrum.push_back(it_mz);
        }
      }
      exp_centroided_.addSpectrum(std::move(spectrum));
    }
    exp_centroided_.updateRanges();
    exp_centroided_.sortSpectra();
    
    // initialise blacklist <blacklist_> and its owners
    blacklist_.reserve(exp_centroided_.getNrSpectra());
    blacklist_owner_.reserve(exp_centroided_.getNrSpectra());
    // loop over spectra
    for (const auto &it_rt : exp_centroided_)
    {
      std::vector<int> blacklist_spectrum(it_rt.size(), -1);
      blacklist_.push_back(blacklist_spectrum);
      blacklist_owner_.push_back(std::move(blacklist_spectrum));
    }
    
  }
  
  MSExperiment& MultiplexFiltering::getCentroidedExperiment()
  {
    return exp_centroided_;
  }

  void MultiplexFiltering::updateWhiteMSExperiment_()
  {
    // reset both the white MS experiment and the corresponding mapping to the complete i.e. original MS experiment
    exp_centroided_white_.clear(true);
    const SignedSize n_spectra = static_cast<SignedSize>(exp_centroided_.size());
    std::vector<MSSpectrum> spectra(n_spectra);
    exp_centroided_mapping_.assign(n_spectra, {});

    // loop over spectra (independent of each other)
    #pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize idx_rt = 0; idx_rt < n_spectra; ++idx_rt)
    {
      const MSSpectrum& spectrum = exp_centroided_[idx_rt];
      const std::vector<int>& blacklist_spectrum = blacklist_[idx_rt];
      MSSpectrum& spectrum_picked_white = spectra[idx_rt];
      std::vector<int>& mapping_spectrum = exp_centroided_mapping_[idx_rt];
      spectrum_picked_white.setRT(spectrum.getRT());

      const Size n_white = static_cast<Size>(std::count(blacklist_spectrum.begin(), blacklist_spectrum.end(), -1));
      spectrum_picked_white.reserve(n_white);
      mapping_spectrum.reserve(n_white);
      // loop over m/z
      for (Size idx_mz = 0; idx_mz < spectrum.size(); ++idx_mz)
      {
        if (blacklist_spectrum[idx_mz] == -1)
        {
          spectrum_picked_white.push_back(spectrum[idx_mz]);
          mapping_spectrum.push_back(static_cast<int>(idx_mz));
        }
      }
    }
    exp_centroided_white_.reserve(n_spectra);
    for (MSSpectrum& spectrum : spectra)
    {
      exp_centroided_white_.addSpectrum(std::move(spectrum));
    }
    exp_centroided_white_.updateRanges();
  }
  
  int MultiplexFiltering::checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak,
                                                   double fraction) const
  {
    // Check that there is a peak.
    int mz_idx = it_rt->findNearest(mz, mz_tolerance);
    if (mz_idx != -1)
    {
      // Any peak with an intensity greater than <fraction>*<intensity_first_peak> is significant.
      MSSpectrum::ConstIterator it_mz = it_rt->begin();
      std::advance(it_mz, mz_idx);
      double intensity = it_mz->getIntensity();

      // Check that the peak is significant.
      if (intensity > fraction * intensity_first_peak)
      {
        // There is a high-intensity peak at the position mz.
        return mz_idx;
      }
    }

    return -1;
  }
  
  bool MultiplexFiltering::filterPeakPositions_(double mz, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak, std::vector<BlacklistEntry>* blacklist_reads) const
  {    
    // check if peak position is blacklisted
    // i.e. -1 = white or 0 = mono-isotopic peak of the lightest (or only) peptide are ok.
    if (blacklist_reads) { blacklist_reads->push_back({peak.getRTidx(), peak.getMZidx(), blacklist_[peak.getRTidx()][peak.getMZidx()]}); }
    if (blacklist_[peak.getRTidx()][peak.getMZidx()] > 0)
    {
      return false;
    }
    
    // determine absolute m/z tolerance in Th
    double mz_tolerance;
    if (mz_tolerance_unit_in_ppm_)
    {
      // m/z tolerance in ppm
      // Note that the absolute tolerance varies minimally within an m/z pattern.
      // Hence we calculate it only once here.
      mz_tolerance = mz * mz_tolerance_ * 1e-6;
    }
    else
    {
      // m/z tolerance in Th
      mz_tolerance = mz_tolerance_;
    }

    // The mass traces of the peptide(s) form a m/z shift pattern. Starting with the mono-isotopic mass trace of each peptide,
    // how long is the series of m/z shifts until the first expected mass trace is missing? We want to see
    // at least isotopes_per_peptide_min_ of these m/z shifts in each peptide. Note that we need to demand subsequent(!) mass traces
    // to be present. Otherwise it would be easy to mistake say a 2+ peptide for a 4+ peptide.
    size_t length = 0;
    bool interrupted = false;
    // The satellites found (spectrum index, peak index, m/z shift index), in this order added to the peak once the
    // pattern is long enough. (Most peaks fail before, and adding to the satellites of a peak allocates.)
    static thread_local std::vector<std::tuple<size_t, size_t, size_t>> found_satellites_of_thread;
    std::vector<std::tuple<size_t, size_t, size_t>>& found_satellites = found_satellites_of_thread;
    found_satellites.clear();
    auto isSatellite = [&peak, &found_satellites](size_t rt_idx, size_t mz_idx)
    {
      return peak.checkSatellite(rt_idx, mz_idx)
             || std::any_of(found_satellites.begin(), found_satellites.end(),
                            [&](const auto& satellite) { return std::get<0>(satellite) == rt_idx && std::get<1>(satellite) == mz_idx; });
    };
    // where the last search of each m/z shift in each spectrum of the band ended on this thread (see findNearestFrom())
    static thread_local std::vector<Size> cursors_of_thread;
    std::vector<Size>& cursors = cursors_of_thread;
    const Size n_shifts = pattern.getMassShiftCount() * isotopes_per_peptide_max_; // (the m/z shift indices used below)
    if (cursors.size() < static_cast<Size>(it_rt_band_end - it_rt_band_begin) * n_shifts)
    {
      cursors.resize(static_cast<Size>(it_rt_band_end - it_rt_band_begin) * n_shifts, 0);
    }
    // loop over isotopes i.e. mass traces within the peptide
    for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
    {
      bool found_in_all_peptides = true;
      
      // loop over peptides
      for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
      {
        // calculate m/z shift index in pattern
        size_t mz_shift_idx = peptide * isotopes_per_peptide_max_ + isotope;
        double mz_shift = pattern.getMZShiftAt(mz_shift_idx);

        bool found = false;
        // loop over spectra in RT band
        for (MSExperiment::ConstIterator it_rt = it_rt_band_begin; it_rt < it_rt_band_end; ++it_rt)
        {
          int i = findNearestFrom(*it_rt, mz + mz_shift, mz_tolerance, cursors[static_cast<Size>(it_rt - it_rt_band_begin) * n_shifts + mz_shift_idx]);
         
          if (i != -1)
          {
            // Note that as primary peaks, satellite peaks are also restricted by the blacklist.
            // The peak can either be pure white i.e. untouched, or have been seen earlier as part of the same mass trace.
            size_t rt_idx = it_rt - it_rt_begin;
            size_t mz_idx = exp_centroided_mapping_.at(it_rt - it_rt_begin).at(i);
            
            // Check that the peak has not been blacklisted and is not already in the satellite set.
            if (blacklist_reads) { blacklist_reads->push_back({rt_idx, mz_idx, blacklist_[rt_idx][mz_idx]}); }
            if (((blacklist_[rt_idx][mz_idx] == -1) || (blacklist_[rt_idx][mz_idx] == static_cast<int>(mz_shift_idx))) && (!isSatellite(rt_idx, mz_idx)))
            {
              found = true;
              found_satellites.emplace_back(rt_idx, mz_idx, mz_shift_idx);
            }
            
          }
        }
        
        if (!found)
        {
          found_in_all_peptides = false;
        }
        
      }
      
      if (found_in_all_peptides && (!interrupted))
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
    for (const auto& [rt_idx, mz_idx, mz_shift_idx] : found_satellites)
    {
      peak.addSatellite(rt_idx, mz_idx, mz_shift_idx);
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
        MSSpectrum::ConstIterator it_mz_first_peak = it_rt->begin();
        std::advance(it_mz_first_peak, mz_idx_first_peak);
        double intensity_first_peak = it_mz_first_peak->getIntensity();

        double mz;
        
        // Check if there is a zeroth peak. The lower isotopologue of a deuterium label lies there, too (setDeuteriumLabels()).
        mz = peak.getMZ() + 2 * pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_) - pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1);
        int mz_idx = checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak, 0.3 + deuteriumLowerFraction_(pattern, peptide, mz_first_peak));
        if (mz_idx != -1)
        {
          // So there is a significant peak to the left. This is only a problem, if this peak is not part of the pattern which we currently detect.
          // (mz_idx is its index in the white spectrum, the satellites hold indices in exp_centroided_. They differ as soon as
          // an earlier pattern has blacklisted a peak of lower m/z in this spectrum.)
          if (!(peak.checkSatellite(peak.getRTidx(), exp_centroided_mapping_.at(peak.getRTidx()).at(mz_idx))))
          {
            return false;
          }
        }

        // Check mistaken charge state identities
        // We are searching the patterns in the order of the most common occurrence (and not decreasing charge state). 
        // That can lead to mistaken charge state identities. Here we check that this is not the case.

        if (pattern.getCharge() == 2)
        {          
          // Is the 2+ pattern really a 4+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)/2 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/2;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak) != -1)
          {
            return false;
          }
          
          // Is the 2+ pattern really a 6+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)*2/3 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/3;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak) != -1)
          {
            return false;
          }
        }

        if (pattern.getCharge() == 3)
        {
          // Is the 3+ pattern really a 6+ pattern?
          mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)/2 + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/2;
          if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak) != -1)
          {
            return false;
          }
        }

        if (pattern.getCharge() == 1)
        {
          for (int c = 2; c < 7; ++c)
          {
            // Is the 1+ pattern really a c+ pattern?
            // (In theory, any charge state c >= 2+ could be mistaken as a 1+. For the sake of run time performance, we only check up to 7+.
            // If we see in any dataset significant number of mistakes for c >= 8+, we will modify this part of the code.)
            mz = peak.getMZ() + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_)*(c-1)/c + pattern.getMZShiftAt(peptide * isotopes_per_peptide_max_ + 1)/c;
            if (checkForSignificantPeak_(mz, 2 * mz_tolerance, it_rt, intensity_first_peak) != -1)
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
  
  void MultiplexFiltering::filterSpectrum_(size_t idx_rt, unsigned pattern_idx, const PeakFilter& filter_peak, MultiplexFilteredMSExperiment& result)
  {
    const SignedSize n_peaks = static_cast<SignedSize>(exp_centroided_white_[idx_rt].size());

#ifdef _OPENMP
    const SignedSize n_threads = omp_get_max_threads();
#else
    const SignedSize n_threads = 1;
#endif
    // A spectrum of few peaks is filtered faster on one thread than the threads take to start. (The serial loop gives the
    // same result.)
    const bool parallel = n_threads > 1 && n_peaks >= std::max<SignedSize>(64, 2 * n_threads);
    if (!parallel)
    {
      for (SignedSize s = 0; s < n_peaks; ++s)
      {
        if (std::optional<MultiplexFilteredPeak> peak = filter_peak(s, nullptr))
        {
          // (blacklisted before it is moved into the result, as the peak it will be there)
          blacklistPeak_(*peak, pattern_idx, result.size());
          result.addPeak(std::move(*peak));
        }
      }
      return;
    }

    // Filter every peak against the blacklist as it is now, and find the entries each one that passes would
    // blacklist. Nothing writes to the blacklist meanwhile.
    // (The buffers are kept from spectrum to spectrum: allocating them anew for every spectrum costs more than filtering
    // its peaks, and every element used is assigned below.)
    std::vector<std::optional<MultiplexFilteredPeak>>& passed = passed_;
    if (blacklist_reads_.size() < static_cast<Size>(n_peaks))
    {
      blacklist_reads_.resize(n_peaks);
      blacklist_writes_.resize(n_peaks);
    }
    if (passed.size() < static_cast<Size>(n_peaks))
    {
      passed.resize(n_peaks);
    }
    // (dynamic: the peaks that pass the first checks take much longer than the others; about 16 chunks per thread)
    const SignedSize chunk = std::max<SignedSize>(1, n_peaks / (16 * n_threads));
    #pragma omp parallel for schedule(dynamic, chunk)
    for (SignedSize s = 0; s < n_peaks; ++s)
    {
      blacklist_reads_[s].clear();
      blacklist_writes_[s].clear();
      passed[s] = filter_peak(s, &blacklist_reads_[s]);
      if (passed[s]) { blacklistEntries_(*passed[s], pattern_idx, blacklist_writes_[s]); }
    }

    // Accept the peaks in order of m/z. A peak is filtered again if a peak accepted before it in this spectrum has
    // changed an entry its decision read, so that every decision sees the blacklist a serial loop would.
    bool blacklisted = false;
    for (SignedSize s = 0; s < n_peaks; ++s)
    {
      const std::vector<BlacklistEntry>& reads = blacklist_reads_[s];
      if (blacklisted
          && std::any_of(reads.begin(), reads.end(),
                         [this](const BlacklistEntry& read) { return blacklist_[read.rt_idx][read.mz_idx] != read.value; }))
      {
        passed[s] = filter_peak(s, nullptr);
        blacklist_writes_[s].clear();
        if (passed[s]) { blacklistEntries_(*passed[s], pattern_idx, blacklist_writes_[s]); }
      }
      if (passed[s])
      {
        result.addPeak(std::move(*passed[s]));
        applyBlacklistEntries_(blacklist_writes_[s], pattern_idx, result.size() - 1);
        blacklisted = true;
      }
    }
  }

  void MultiplexFiltering::blacklistPeak_(const MultiplexFilteredPeak& peak, unsigned pattern_idx, size_t peak_idx)
  {
    std::vector<BlacklistEntry> entries;
    blacklistEntries_(peak, pattern_idx, entries);
    applyBlacklistEntries_(entries, pattern_idx, peak_idx);
  }

  void MultiplexFiltering::applyBlacklistEntries_(const std::vector<BlacklistEntry>& entries, unsigned pattern_idx, size_t peak_idx)
  {
    // owners and records are indexed by int (see blacklist_owner_); each entry adds at most one record
    constexpr Size max_count = static_cast<Size>(std::numeric_limits<int>::max());
    if (filtered_peaks_.size() >= max_count || peak_idx >= max_count || entries.size() > max_count - blacklist_owner_records_.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filtered_peaks_.size() + 1,
                                   "Too many filtered peaks, or peaks reserved by more than one of them, for the blacklist owners (limit: INT_MAX each).");
    }
    const int owner = static_cast<int>(filtered_peaks_.size());
    filtered_peaks_.emplace_back(pattern_idx, static_cast<unsigned>(peak_idx));
    for (const BlacklistEntry& entry : entries)
    {
      int& value = blacklist_[entry.rt_idx][entry.mz_idx];
      int& owned_by = blacklist_owner_[entry.rt_idx][entry.mz_idx];
      if (owned_by == -1)
      {
        owned_by = owner;
      }
      else if (owned_by >= 0)
      {
        // a peak another filtered peak reserved before keeps a record of the first one, and of the last one from now on
        if (owned_by != owner)
        {
          blacklist_owner_records_.push_back({owned_by, value, owner});
          owned_by = -2 - static_cast<int>(blacklist_owner_records_.size() - 1);
        }
      }
      else
      {
        blacklist_owner_records_[static_cast<Size>(-2 - owned_by)].last = owner;
      }
      value = entry.value;
    }
  }

  void MultiplexFiltering::blacklistEntries_(const MultiplexFilteredPeak& peak, unsigned pattern_idx, std::vector<BlacklistEntry>& entries) const
  {
    // determine absolute m/z tolerance in Th
    double mz_tolerance;
    if (mz_tolerance_unit_in_ppm_)
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
    
    // Determine the RT boundaries for each of the mass traces.
    const std::multimap<size_t, MultiplexSatelliteCentroided >& satellites = peak.getSatellites();
    // <rt_boundaries> is a map from the mass trace index to the spectrum indices for beginning and end of the mass trace.
    std::map<size_t, std::pair<size_t, size_t> > rt_boundaries;
    // loop over satellites
    for (const auto &it : satellites)
    {
      size_t idx_masstrace = it.first;    // mass trace index i.e. the index within the peptide multiplet pattern
      if (!rt_boundaries.contains(idx_masstrace))
      {
        // That's the first satellite within this mass trace.
        rt_boundaries[idx_masstrace] = std::make_pair((it.second).getRTidx(), (it.second).getRTidx());
      }
      else
      {
        // We have seen a satellite of this mass trace before.
        size_t idx_min = std::min((it.second).getRTidx(), rt_boundaries[idx_masstrace].first);
        size_t idx_max = std::max((it.second).getRTidx(), rt_boundaries[idx_masstrace].second);
        
        rt_boundaries[idx_masstrace] = std::make_pair(idx_min, idx_max);
      }
    }
    
    // where the last search of each mass trace in each spectrum ended on this thread (see findNearestFrom(); the peaks
    // of a spectrum are blacklisted in order of m/z)
    static thread_local std::vector<Size> cursors_of_thread;
    std::vector<Size>& cursors = cursors_of_thread;
    const Size n_traces = patterns_[pattern_idx].getMassShiftCount() * isotopes_per_peptide_max_; // (the mass trace indices)
    if (cursors.size() < exp_centroided_.size() * n_traces)
    {
      cursors.resize(exp_centroided_.size() * n_traces, 0);
    }

    // Blacklist all peaks along the mass traces
    // loop over mass traces (i.e. the mass trace boundaries)
    for (const auto &it : rt_boundaries)
    {
      double mz = peak.getMZ() + patterns_[pattern_idx].getMZShiftAt(it.first);
      
      // Extend the RT boundary by rt_band_ earlier
      MSExperiment::ConstIterator it_rt_begin = exp_centroided_.begin() + (it.second).first;
      it_rt_begin = exp_centroided_.RTBegin(it_rt_begin->getRT() - 2 * rt_band_);
      
      // Extend the RT boundary by rt_band_ later
      MSExperiment::ConstIterator it_rt_end = exp_centroided_.begin() + (it.second).second;
      it_rt_end = exp_centroided_.RTBegin(it_rt_end->getRT() + 2 * rt_band_);
      
      // prepare for loop
      if (it_rt_end != exp_centroided_.end())
      {
        ++it_rt_end;
      }
      
      // loop over RT along the mass trace
      for (MSExperiment::ConstIterator it_rt = it_rt_begin; it_rt < it_rt_end; ++it_rt)
      {
        int idx_mz = findNearestFrom(*it_rt, mz, mz_tolerance, cursors[static_cast<Size>(it_rt - exp_centroided_.begin()) * n_traces + it.first]);
        
        if (idx_mz != -1)
        {
          // blacklist entries: -1 = white, any isotope pattern index (it.first) = black
          entries.push_back({static_cast<size_t>(it_rt - exp_centroided_.begin()), static_cast<size_t>(idx_mz), static_cast<int>(it.first)});
        }
      }
      
    }
    
  }
  
  MSExperiment MultiplexFiltering::getBlacklist()
  {
    MSExperiment exp_blacklist;
    
    // loop over spectra
    for (const auto &it_rt : exp_centroided_)
    {
      MSSpectrum spectrum_black;
      spectrum_black.setRT(it_rt.getRT());
      
      // loop over m/z
      for (const auto &it_mz : it_rt)
      {
        // transfer all peaks which are not white (i.e. not -1)
        if (blacklist_[&it_rt - &exp_centroided_[0]][&it_mz - &it_rt[0]] != -1)
        {
          spectrum_black.push_back(it_mz);
        }
      }
      exp_blacklist.addSpectrum(spectrum_black);
    }
    exp_blacklist.updateRanges();
    
    return exp_blacklist;
  }

  std::vector<std::vector<MultiplexFiltering::BlacklistOwner>> MultiplexFiltering::getBlacklistOwners() const
  {
    std::vector<std::vector<BlacklistOwner>> owners;
    owners.reserve(blacklist_.size());
    for (Size rt_idx = 0; rt_idx < blacklist_.size(); ++rt_idx)
    {
      owners.push_back(getBlacklistOwners(rt_idx));
    }
    return owners;
  }

  std::vector<MultiplexFiltering::BlacklistOwner> MultiplexFiltering::getBlacklistOwners(Size spectrum_index) const
  {
    if (spectrum_index >= blacklist_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, static_cast<SignedSize>(spectrum_index), blacklist_.size());
    }
    auto writer = [this](int owner, int value) -> BlacklistWriter
    {
      if (owner < 0) { return {std::numeric_limits<unsigned>::max(), std::numeric_limits<unsigned>::max(), value}; }
      return {filtered_peaks_[owner].first, filtered_peaks_[owner].second, value};
    };
    const std::vector<int>& values = blacklist_[spectrum_index];
    const std::vector<int>& owners_of_peaks = blacklist_owner_[spectrum_index];
    // the same peaks in the same order as getBlacklist(), which keeps the peaks that are not white
    std::vector<BlacklistOwner> owners;
    owners.reserve(values.size() - static_cast<Size>(std::count(values.begin(), values.end(), -1)));
    for (Size mz_idx = 0; mz_idx < values.size(); ++mz_idx)
    {
      const int value = values[mz_idx];
      if (value == -1) { continue; }
      const int owner = owners_of_peaks[mz_idx];
      if (owner <= -2)
      {
        const BlacklistOwnerRecord& record = blacklist_owner_records_[static_cast<Size>(-2 - owner)];
        owners.push_back({writer(record.first, record.first_value), writer(record.last, value)});
      }
      else
      {
        const BlacklistWriter last = writer(owner, value);
        owners.push_back({last, last});
      }
    }
    return owners;
  }
  
  const IsotopeDistribution& MultiplexFiltering::averagineDistribution_(double mass) const
  {
    if (averagine_type_ == "peptide")
    {
      // the same as CoarseIsotopePatternGenerator(isotopes_per_peptide_max_).estimateFromPeptideWeight(mass), computed
      // once per element composition of the averagine model (many masses share one)
      return DIAHelpers::averagineIsotopeDistribution(mass, static_cast<int>(isotopes_per_peptide_max_));
    }
    static thread_local IsotopeDistribution distribution;
    CoarseIsotopePatternGenerator solver(isotopes_per_peptide_max_);
    if (averagine_type_ == "RNA")
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
    return distribution;
  }

  void MultiplexFiltering::setDeuteriumLabels(const std::map<std::string, unsigned>& deuterium_per_label, double enrichment)
  {
    if (!(enrichment > 0.0 && enrichment <= 1.0))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The 2H enrichment of deuterium labels must lie in (0, 1].",
                                    std::to_string(enrichment));
    }
    deuterium_enrichment_ = enrichment;
    deuterium_.clear();
    for (const MultiplexIsotopicPeakPattern& pattern : patterns_)
    {
      const int index = pattern.getMassShiftIndex();
      if (index < 0) { continue; }
      std::vector<unsigned> deuterium;
      bool any = false;
      const MultiplexDeltaMasses mass_shifts = pattern.getMassShifts(); // a copy: getDeltaMasses() refers into it
      for (const MultiplexDeltaMasses::DeltaMass& delta_mass : mass_shifts.getDeltaMasses())
      {
        unsigned n = 0;
        for (const std::string& label : delta_mass.label_set)
        {
          const auto it = deuterium_per_label.find(label);
          if (it != deuterium_per_label.end()) { n += it->second; }
        }
        deuterium.push_back(n);
        any = any || n > 0;
      }
      if (!any) { continue; }
      if (deuterium_.size() <= static_cast<size_t>(index)) { deuterium_.resize(index + 1); }
      deuterium_[index] = std::move(deuterium);
    }
  }

  double MultiplexFiltering::deuteriumLowerFraction(unsigned deuterium, double enrichment, const IsotopeDistribution& isotopes)
  {
    if (deuterium == 0 || !(enrichment > 0.0 && enrichment < 1.0) || isotopes.size() < 2) { return 0.0; }
    const double odds = (1.0 - enrichment) / enrichment;
    double q = 1.0; // q_j, relative to the fully labelled molecules (j = 0)
    double lower = 0.0, monoisotopic = isotopes[0].getIntensity();
    for (unsigned j = 1; j <= deuterium && j <= isotopes.size(); ++j)
    {
      q *= odds * static_cast<double>(deuterium - j + 1) / static_cast<double>(j);
      lower += q * isotopes[j - 1].getIntensity();
      if (j < isotopes.size()) { monoisotopic += q * isotopes[j].getIntensity(); }
    }
    return monoisotopic > 0.0 ? lower / monoisotopic : 0.0;
  }

  double MultiplexFiltering::deuteriumLowerFraction_(const MultiplexIsotopicPeakPattern& pattern, size_t peptide, double mz) const
  {
    const int index = pattern.getMassShiftIndex();
    if (index < 0 || static_cast<size_t>(index) >= deuterium_.size() || peptide >= deuterium_[index].size() || deuterium_[index][peptide] == 0)
    {
      return 0.0;
    }
    return deuteriumLowerFraction(deuterium_[index][peptide], deuterium_enrichment_, averagineDistribution_(mz * pattern.getCharge()));
  }

  bool MultiplexFiltering::filterAveragineModel_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const
  {
    // construct averagine distribution
    double mass = peak.getMZ() * pattern.getCharge();
    const IsotopeDistribution& distribution = averagineDistribution_(mass);
    
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
        std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites;
        satellites = peak.getSatellites().equal_range(idx);
              
        int count = 0;
        double sum_intensities = 0;
        
        // loop over satellites in mass trace
        for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = satellites.first; satellite_it != satellites.second; ++satellite_it)
        {
          // find indices of the peak
          size_t rt_idx = (satellite_it->second).getRTidx();
          size_t mz_idx = (satellite_it->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt = exp_centroided_.begin();
          std::advance(it_rt, rt_idx);
          MSSpectrum::ConstIterator it_mz = it_rt->begin();
          std::advance(it_mz, mz_idx);
                    
          ++count;
          sum_intensities += it_mz->getIntensity();
        }
        
        if (count > 0)
        {
          //intensities_model.push_back(distribution.getContainer()[isotope].second);
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
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0, "intensity vectors too small for correlation calculation");
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

  bool MultiplexFiltering::filterPeptideCorrelation_(const MultiplexIsotopicPeakPattern& pattern, const MultiplexFilteredPeak& peak) const
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
                    
          std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_1;
          std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_2;
          satellites_1 = peak.getSatellites().equal_range(idx_1);
          satellites_2 = peak.getSatellites().equal_range(idx_2);
          
          // loop over satellites in mass trace 1
          for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_1 = satellites_1.first; satellite_it_1 != satellites_1.second; ++satellite_it_1)
          {
            size_t rt_idx_1 = (satellite_it_1->second).getRTidx();
  
            // loop over satellites in mass trace 2
            for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_2 = satellites_2.first; satellite_it_2 != satellites_2.second; ++satellite_it_2)
            {
              size_t rt_idx_2 = (satellite_it_2->second).getRTidx();

              if (rt_idx_1 == rt_idx_2)
              {
                size_t mz_idx_1 = (satellite_it_1->second).getMZidx();
                size_t mz_idx_2 = (satellite_it_2->second).getMZidx();
                
                // find peak itself
                MSExperiment::ConstIterator it_rt_1 = exp_centroided_.begin();
                MSExperiment::ConstIterator it_rt_2 = exp_centroided_.begin();
                std::advance(it_rt_1, rt_idx_1);
                std::advance(it_rt_2, rt_idx_2);
                MSSpectrum::ConstIterator it_mz_1 = it_rt_1->begin();
                MSSpectrum::ConstIterator it_mz_2 = it_rt_2->begin();
                std::advance(it_mz_1, mz_idx_1);
                std::advance(it_mz_2, mz_idx_2);
  
                intensities_1.push_back(it_mz_1->getIntensity());
                intensities_2.push_back(it_mz_2->getIntensity());
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
