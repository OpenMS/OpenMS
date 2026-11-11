// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include "ProSEAlgorithmImpl.h"

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <thread>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  namespace
  {
    // The peak list and native ID of @p from into the empty spectrum @p to (and, if @p with_precursors, its MS level and
    // precursors): a full MSSpectrum copy also copies the meta values, the data arrays and the other settings, heap blocks
    // per spectrum that the readers of these copies never look at (the native ID identifies the copy).
    void copyPeakList_(const MSSpectrum& from, MSSpectrum& to, bool with_precursors)
    {
      to.reserve(from.size());
      to.insert(to.end(), from.begin(), from.end());
      to.setNativeID(from.getNativeID());
      if (with_precursors)
      {
        to.setMSLevel(from.getMSLevel());
        to.setPrecursors(from.getPrecursors());
      }
    }
  }

  ProSEAlgorithm::SpectrumGenerators_ ProSEAlgorithm::spectrumGenerators_() const
  {
    // Scoring, annotation and calibration use the ion series of the fragment index, so that
    // a candidate is scored against the ions it was retrieved by.
    auto configure = [this](TheoreticalSpectrumGenerator& tsg, bool electron)
    {
      Param p(tsg.getParameters());
      p.setValue("add_first_prefix_ion", "true");
      p.setValue("add_metainfo", "true");
      p.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
      p.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
      p.setValue("add_c_ions", (add_c_ions_ || electron) ? "true" : "false");
      p.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
      p.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
      p.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
      p.setValue("add_zp1_ions", (add_zp1_ions_ || electron) ? "true" : "false");
      tsg.setParameters(p);
    };
    SpectrumGenerators_ generators;
    configure(generators.standard, false);
    configure(generators.electron, true);
    generators.by_activation = ions_by_activation_;
    return generators;
  }

  std::vector<double> ProSEAlgorithm::localPeakDensities_(const MSSpectrum& spectrum)
  {
    // Both bounds only move forward: O(number of peaks), once per spectrum.
    constexpr double half_window = 50.0;
    std::vector<double> densities(spectrum.size());
    Size left = 0, right = 0;
    for (Size i = 0; i < spectrum.size(); ++i)
    {
      const double mz = spectrum[i].getMZ();
      while (left < spectrum.size() && spectrum[left].getMZ() < mz - half_window) ++left;
      while (right < spectrum.size() && spectrum[right].getMZ() <= mz + half_window) ++right;
      densities[i] = static_cast<double>(right - left) / (2.0 * half_window);
    }
    return densities;
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical,
      const std::vector<double>& densities, double tolerance, bool ppm)
  {
    std::vector<double> alternatives;
    return localFragmentEvidence_(spectrum, theoretical, std::numeric_limits<int>::max(), densities, tolerance, ppm, alternatives);
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical, int max_charge,
      const std::vector<double>& densities, double tolerance, bool ppm, std::vector<double>& alternatives)
  {
    if (spectrum.empty() || theoretical.empty()) return LocalFragmentEvidence_();
    OPENMS_PRECONDITION(!theoretical.getIntegerDataArrays().empty()
      && theoretical.getIntegerDataArrays()[0].size() == theoretical.size(), "Aligned fragment charges required")
    return localFragmentEvidence_(spectrum, theoretical, theoretical.getIntegerDataArrays()[0], max_charge, densities, tolerance,
                                  ppm, alternatives);
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical, const std::vector<Int>& charges, int max_charge,
      const std::vector<double>& densities, double tolerance, bool ppm, std::vector<double>& alternatives)
  {
    // The ions of charge <= max_charge of a TSG spectrum of charges 1..Z are the TSG spectrum of charges
    // 1..max_charge, in the same order: the generator adds each charge independently and sorts stably.
    // Annotation therefore passes its own theoretical spectrum instead of generating a second one.
    LocalFragmentEvidence_ evidence;
    if (spectrum.empty() || theoretical.empty()) return evidence;
    OPENMS_PRECONDITION(densities.size() == spectrum.size(), "One local density per experimental peak required")
    OPENMS_PRECONDITION(charges.size() == theoretical.size(), "Aligned fragment charges required")

    static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
    static const double ammonia_mass = EmpiricalFormula("NH3").getMonoWeight();
    alternatives.clear();
    alternatives.reserve(3 * theoretical.size());
    for (Size i = 0; i < theoretical.size(); ++i)
    {
      OPENMS_PRECONDITION(charges[i] > 0, "Positive fragment charge required")
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      alternatives.push_back(mz);
      for (const double loss : {water_mass, ammonia_mass})
      {
        const double loss_mz = mz - loss / charges[i];
        if (loss_mz > 0.0) alternatives.push_back(loss_mz);
      }
    }
    std::sort(alternatives.begin(), alternatives.end());

    for (Size i = 0; i < theoretical.size(); ++i)
    {
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      const double width = ppm ? mz * tolerance * 1e-6 : tolerance;
      const Size nearest = spectrum.findNearest(mz);
      const double observed = spectrum[nearest].getMZ();
      if (std::abs(mz - observed) > width) continue;

      // The matched peak itself is in its density window, hence density > 0.
      const double density = densities[nearest];
      evidence.chance_match_surprise += std::max(0.0, -std::log(2.0 * width * density));
      const auto first = std::lower_bound(alternatives.begin(), alternatives.end(), observed - width);
      const auto last = std::upper_bound(first, alternatives.end(), observed + width);
      // One intact hypothesis is the ion currently being scored. All other
      // hypotheses, including coincident ions, represent alternative assignments.
      const Size count = static_cast<Size>(last - first);
      const Size competing = count > 0 ? count - 1 : 0;
      evidence.mass_competition_evidence += 1.0 / (1.0 + competing + density);
    }
    return evidence;
  }

  bool ProSEAlgorithm::isElectronActivated_(const MSSpectrum& spectrum)
  {
    if (spectrum.getPrecursors().empty()) return false;
    for (const Precursor::ActivationMethod method : spectrum.getPrecursors()[0].getActivationMethods())
    {
      if (method == Precursor::ActivationMethod::ETD || method == Precursor::ActivationMethod::ECD
          || method == Precursor::ActivationMethod::EThcD || method == Precursor::ActivationMethod::ETciD)
      {
        return true;
      }
    }
    return false;
  }

  Size ProSEAlgorithm::countElectronActivated_(const PeakMap& spectra) const
  {
    if (!ions_by_activation_) return 0;
    return static_cast<Size>(std::count_if(spectra.begin(), spectra.end(),
                                           [](const MSSpectrum& spectrum) { return isElectronActivated_(spectrum); }));
  }

  namespace
  {
    // Keeps the peaks WindowMower::filterPeakSpectrumForTopNInJumpingWindow() keeps: the peak_count most intense peaks of
    // each m/z window of window_size that starts at a peak (of the last window a share of peak_count by its width), and
    // every peak equal to a kept one. It selects them through indices instead of copies of the peaks, and finds the
    // equal peaks among their neighbours instead of searching all kept peaks for each peak. The same std::partial_sort
    // calls on the same intensities select the same peaks.
    void filterTopNInJumpingWindows(MSSpectrum& spectrum, double window_size, UInt peak_count)
    {
      if (spectrum.empty()) { return; }
      spectrum.sortByPosition();
      const Size n = spectrum.size();
      std::vector<char> kept(n, 0);
      std::vector<Size> window;
      const auto more_intense = [&spectrum](Size a, Size b) { return spectrum[b].getIntensity() < spectrum[a].getIntensity(); };
      const auto keep_most_intense = [&](Size begin, Size end, Size count)
      {
        if (end - begin > count)
        {
          window.resize(end - begin);
          std::iota(window.begin(), window.end(), begin);
          std::partial_sort(window.begin(), window.begin() + count, window.end(), more_intense);
          for (Size k = 0; k < count; ++k) { kept[window[k]] = 1; }
        }
        else
        {
          std::fill(kept.begin() + begin, kept.begin() + end, 1);
        }
      };
      double window_start = spectrum[0].getMZ();
      Size begin = 0;
      for (Size i = 0; i != n; ++i)
      {
        if (spectrum[i].getMZ() - window_start < window_size) { continue; }
        // a gap may leave windows empty: the next window starts at the next peak
        window_start = spectrum[i].getMZ();
        keep_most_intense(begin, i, peak_count);
        begin = i;
      }
      const double last_window_fraction = (spectrum[n - 1].getMZ() - window_start) / window_size;
      keep_most_intense(begin, n, static_cast<Size>(std::round(last_window_fraction * peak_count)));

      // peaks of equal m/z are neighbours (the spectrum is sorted by m/z)
      std::vector<Size> selected;
      selected.reserve(n);
      for (Size i = 0; i < n;)
      {
        Size j = i + 1;
        while (j < n && spectrum[j].getMZ() == spectrum[i].getMZ()) { ++j; }
        for (Size a = i; a < j; ++a)
        {
          for (Size b = i; b < j; ++b)
          {
            if (kept[b] && spectrum[b] == spectrum[a])
            {
              selected.push_back(a);
              break;
            }
          }
        }
        i = j;
      }
      spectrum.select(selected);
    }
  }

  void ProSEAlgorithm::filterLocalPeaks_(MSSpectrum& spectrum, Size peaks_per_window)
  {
    // Work on indices so all peak-associated data arrays survive the selection.
    // An isolated high-m/z ion is still evidence; its window's observed width
    // must not round the retention quota down to zero.
    std::vector<Size> indices(spectrum.size()), selected;
    std::iota(indices.begin(), indices.end(), 0);
    selected.reserve(spectrum.size());
    for (Size begin = 0; begin < spectrum.size();)
    {
      Size end = begin + 1;
      while (end < spectrum.size() && spectrum[end].getMZ() - spectrum[begin].getMZ() < 100.0)
      {
        ++end;
      }
      const Size keep = std::min(peaks_per_window, end - begin);
      std::partial_sort(indices.begin() + begin, indices.begin() + begin + keep, indices.begin() + end, [&spectrum](Size a, Size b) {
        if (spectrum[a].getIntensity() != spectrum[b].getIntensity()) { return spectrum[a].getIntensity() > spectrum[b].getIntensity(); }
        return a < b;
      });
      selected.insert(selected.end(), indices.begin() + begin, indices.begin() + begin + keep);
      begin = end;
    }
    std::sort(selected.begin(), selected.end());
    spectrum.select(selected);
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type)
  {
    preprocessSpectra_(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, deisotope_requested, peaks_keep_n,
                       peaks_window_top, window_type, DeisotopingSettings_{});
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type,
                                          const DeisotopingSettings_& deisotoping,
                                          FragmentIonLikelihoodModel::PeakLists* ion_evidence,
                                          bool ion_evidence_scored_peaks,
                                          PeakMap* evidence_spectra,
                                          PeakMap* query_spectra)
  {
    // Intensity threshold + normalization used to run here as two extra SERIAL full-map
    // passes. Both are strictly per-spectrum: ThresholdMower::filterPeakMap and
    // Normalizer::filterPeakMap are literally "for (auto& s : exp) filterSpectrum(s);" and
    // neither iterates chromatograms. They are therefore applied at the top of the parallel
    // loop below instead, which is per-spectrum equivalent and removes two full sweeps over
    // the peak data. Both objects are configured once here; like nlargest_filter below, each
    // OpenMP thread works on its own copy (made before the loop): ThresholdMower stores its 'threshold'
    // Param in a member on every call, and concurrent writes are a data race even when every
    // thread writes the same value. One copy per thread costs a few Param copies per search.
    // Peaks without intensity (zero or negative, e.g. empty centroids) would still count as
    // matched ions, so they are removed. Nothing else is (every positive float intensity,
    // subnormal ones included, passes): the ThresholdMower default of 0.05 is an absolute
    // cutoff on the raw intensities before normalization and deleted real peaks from
    // intensity-scaled input (e.g. spectra normalized to a base peak of 1).
    ThresholdMower threshold_mower_filter;
    Param threshold_param = threshold_mower_filter.getParameters();
    threshold_param.setValue("threshold", static_cast<double>(std::numeric_limits<float>::denorm_min()));
    threshold_mower_filter.setParameters(threshold_param);
    Normalizer normalizer;

    // sort by rt; done before the loop because sortSpectra(false) only permutes whole
    // spectra (it does not touch peak order, see MSExperiment::sortSpectra) and therefore
    // commutes with the per-spectrum filters that now run inside the loop.
    if (!inRetentionTimeOrder(exp)) { exp.sortSpectra(false); } // the spectra of loadMS2Spectra() are sorted already
    if (ion_evidence != nullptr) { ion_evidence->reset(exp.size()); }
    if (evidence_spectra != nullptr)
    {
      evidence_spectra->clear(true);
      resizeWithParallelFaults(evidence_spectra->getSpectra(), exp.size());
    }
    if (query_spectra != nullptr)
    {
      query_spectra->clear(true);
      resizeWithParallelFaults(query_spectra->getSpectra(), exp.size());
    }

    // filter settings: the most intense peaks_window_top peaks per 100 Th window (jumping windows as WindowMower's,
    // unless full_window_quota)
    const bool full_window_quota
      = window_type == "jump_full"
        || (window_type == "auto" && Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm));

    // Resolution-aware peak retention. peaks_keep_n == 0 => auto. For HIGH-resolution fragments
    // (within the deisotoper range, <= 0.1 Da / <= 100 ppm) keep the legacy cap of 400. For
    // LOW-resolution data (e.g. 0.5 Da ion-trap CID — the same regime where deisotoping is
    // skipped) a wide match window admits many spurious low-intensity matches that inflate the
    // count-sensitive HyperScore for targets AND decoys alike, collapsing the target/decoy
    // margin; retain far fewer peaks. Random matches scale with peak density x tolerance width,
    // so sqrt-dampen around the 0.02 Da reference and clamp to [60, 400]:
    //   0.2 Da -> 126, 0.5 Da -> 80, 1 Da -> 60.
    Size effective_keep_n = peaks_keep_n;
    if (effective_keep_n == 0)
    {
      if (Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm))
      {
        effective_keep_n = 400; // high-resolution: unchanged legacy behavior
      }
      else
      {
        const double eff_tol_da = fragment_mass_tolerance_unit_ppm
          ? fragment_mass_tolerance * 500.0 * 1e-6   // ppm -> Da at a 500 m/z reference
          : fragment_mass_tolerance;
        const double n = std::round(400.0 * std::sqrt(0.02 / std::max(eff_tol_da, 1e-6)));
        effective_keep_n = static_cast<Size>(std::min(400.0, std::max(60.0, n)));
      }
    }
    NLargest nlargest_filter = NLargest(effective_keep_n);

    // Deisotope only when requested (param fragment:deisotope, resolved in
    // updateMembers_) AND the tolerance is in the Deisotoper's supported range.
    // The range guard keeps the Deisotoper call below from ever throwing, so no
    // exception can escape the OpenMP region (OpenMS#9619). Mode resolution and the
    // low-resolution warning are handled once, in updateMembers_.
    const bool do_deisotope = deisotope_requested &&
      Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);

    // The per-thread copies of the filters are made here: copied on entry to the parallel region (firstprivate), they
    // were outside any exception handling, and an exception from a copy (std::bad_alloc) terminated the process.
#ifdef _OPENMP
    const Size num_threads = static_cast<Size>(omp_get_max_threads());
#else
    const Size num_threads = 1;
#endif
    std::vector<ThresholdMower> threshold_mower_filters(num_threads, threshold_mower_filter);
    std::vector<Normalizer> normalizers(num_threads, normalizer);
    std::vector<NLargest> nlargest_filters(num_threads, nlargest_filter);

    // An exception on a worker thread must not leave the parallel region, which would terminate the process: the
    // first one is rethrown after the loop.
    std::exception_ptr preprocessing_error;
    std::atomic<bool> preprocessing_failed{false};
    // (static: freeInParallel() frees the copies made here in the same parts, each in its own thread's arena)
#pragma omp parallel for schedule(static) default(none) shared(exp, evidence_spectra, query_spectra, do_deisotope, fragment_mass_tolerance, \
                                                fragment_mass_tolerance_unit_ppm, full_window_quota, peaks_window_top, \
                                                deisotoping, ion_evidence, ion_evidence_scored_peaks, \
                                                threshold_mower_filters, normalizers, nlargest_filters, \
                                                preprocessing_error, preprocessing_failed)
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      if (preprocessing_failed.load(std::memory_order_relaxed)) continue;
      try
      {
#ifdef _OPENMP
      const Size thread = static_cast<Size>(omp_get_thread_num());
#else
      const Size thread = 0;
#endif
      // remove 0 intensities, then normalize (formerly two serial full-map passes)
      threshold_mower_filters[thread].filterPeakSpectrum(exp[exp_index]);
      normalizers[thread].filterPeakSpectrum(exp[exp_index]);

      // sort by mz
      exp[exp_index].sortByPosition();

      // fragment:query_spectrum=raw retrieves candidates with this peak list (before deisotoping and local filtering).
      // (only what the candidate query reads: the peaks, MS level and precursors; not the meta data and data arrays)
      if (query_spectra != nullptr) { copyPeakList_(exp[exp_index], (*query_spectra)[exp_index], true); }

      // deisotope (skipped for low-resolution data; see do_deisotope above)
      // Isotope intensities must fall from the monoisotopic peak on (start_intensity_check = 1).
      // With the library default of 2, a small peak one isotope spacing below a fragment ion
      // became the envelope's monoisotopic peak and the ion itself was removed as its isotope.
      // TMT/TMTpro-labelled fragments carry such a peak (reagent isotope impurity), and dense
      // Orbitrap Astral and timsTOF spectra often hold one by chance.
      // The envelope rule is set by fragment:deisotope_min_peaks, _charge_cap and _sum_intensity; their defaults
      // are Sage-like (two peaks, charges up to the precursor charge but at most 3, summed intensity; the earlier
      // rule was three peaks, charges 1-3, own intensity). The Sage-like rule differs from Sage itself: Sage links isotope
      // pairs with a fixed 10 ppm tolerance, requires each isotope peak to be weaker than its parent and does not
      // cap the charge at 3, while the OpenMS deisotoper extends an envelope from its monoisotopic peak with the
      // search tolerance and stops at the first isotope peak that is more intense than its predecessor.
      // The deisotoper does not exclude peaks that already belong to an earlier envelope, so an isotope peak can be
      // claimed by two monoisotopic peaks and, with summing, counted in both (two-peak envelopes make this more
      // frequent; Sage likewise adds a peak to several parents of the same charge). It does not depend on labels.
      if (do_deisotope)
      {
        int max_charge = 3;
        if (deisotoping.charge_cap_precursor && ! exp[exp_index].getPrecursors().empty())
        {
          const int precursor_charge = exp[exp_index].getPrecursors()[0].getCharge();
          if (precursor_charge > 0) { max_charge = std::min(max_charge, precursor_charge); }
        }
        Deisotoper::deisotopeAndSingleCharge(exp[exp_index],
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
          1, max_charge,  // min / max charge
          false,  // keep only deisotoped
          deisotoping.min_peaks, 10,  // min / max isopeaks
          true,   // convert fragment m/z to mono-charge
          false,  // annotate charge
          false,  // annotate isotopic peak counts
          true,   // decreasing isotope intensities
          1,      // start the intensity check at the monoisotopic peak
          deisotoping.sum_intensity);  // the monoisotopic peak carries the envelope's intensity
      }

      // ion priors: every peak after deisotoping, before the window and top-N filters
      if (ion_evidence != nullptr && !ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
      // Local fragment evidence uses this peak list: after deisotoping, so that densities and
      // matching share the m/z space of the scoring spectrum, but before local and top-N filtering.
      // (only what localPeakDensities_ and localFragmentEvidence_ read: the peaks)
      if (evidence_spectra != nullptr) { copyPeakList_(exp[exp_index], (*evidence_spectra)[exp_index], false); }

      // remove noise
      if (full_window_quota) { filterLocalPeaks_(exp[exp_index], static_cast<Size>(peaks_window_top)); }
      else { filterTopNInJumpingWindows(exp[exp_index], 100.0, static_cast<UInt>(peaks_window_top)); }
      nlargest_filters[thread].filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();

      // ion priors on the scored peaks
      if (ion_evidence != nullptr && ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_preprocessing_error)
        {
          if (!preprocessing_error) preprocessing_error = std::current_exception();
        }
        preprocessing_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (preprocessing_error) std::rethrow_exception(preprocessing_error);
  }
} // namespace OpenMS
