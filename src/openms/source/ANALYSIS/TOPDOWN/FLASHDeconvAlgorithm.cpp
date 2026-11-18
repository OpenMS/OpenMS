// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// $Authors: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandJoin.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownIsobaricQuantification.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <tuple>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
inline const Size max_peak_count_for_centroid_ = 3e4;
inline const Size max_peak_count_for_profile_ = 1e5;

FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(): DefaultParamHandler("FLASHDeconvAlgorithm"), ProgressLogger()
{
  // ADR-0046: a commanded MS2 is LOCATED through the acquisition's scan_commands.tsv, not searched for.
  defaults_.setValue("scan_commands", "", "Provide the scan_commands.tsv written by FLASHIda (in the run folder <raw name>_<yyyy-MM-dd-HH-mm-ss>) to couple this mzML with its FLASHIda acquisition. Needs an mzML that keeps each spectrum's 'scan description' (msconvert does). Replaces the former ida_log parameter.");
  defaults_.setValue(
    "report_FDR", "false",
    "Generate q-values (approximate FDR) for deconvolved masses. Decoy masses are also reported which were used for FDR calculation. (Beta version)");
  defaults_.setValidStrings("report_FDR", {"true", "false"});

  defaults_.setValue("allowed_isotope_error", 0,
                     "Tolerance for isotope index errors when calculating FDR. For instance, setting a value of 2 permits the inclusion of up to 2 isotope errors as valid matches. Beta version.");
  defaults_.addTag("allowed_isotope_error", "advanced");

  defaults_.setValue("use_RNA_averagine", "false", "Use the RNA (nucleotide) averagine model for deconvolution.");
  defaults_.setValidStrings("use_RNA_averagine", {"true", "false"});
  defaults_.addTag("use_RNA_averagine", "advanced");

//  defaults_.setValue(
//    "precursor_MS1_window", 1,
//    "Number of MS1 spectra around each MS2 spectrum to search for precursor peaks when determining the MS2 precursors. For MS2 spectrum, the mass of precursor ion should be determined for better deconvolution and reliable identification. "
//    "If the mass of precursor ion is not found in the immediately preceding MS1 spectrum, previous or next MS1 spectra may be used instead. "
//    "This parameter determines up to how many MS1 spectra around each MS2 spectrum will be searched.");
//  defaults_.setMinInt("precursor_MS1_window", 1);
//  defaults_.addTag("precursor_MS1_window", "advanced");

  defaults_.setValue(
    "isolation_window", 5.0,
    "Specify the isolation window width for precursor determination. Used when this information is absent in the mzML file.");
  defaults_.addTag("isolation_window", "advanced");

  defaults_.setValue("precursor_rescue_window", 2,
                     "Number of MS1 spectra before and after the precursor (preceding MS1) spectrum of an MS2 spectrum that are searched for its "
                     "precursor mass when the precursor spectrum has no deconvolved mass in the isolation window: in an LC-MS run the precursor "
                     "elutes over several MS1 scans, so a mass missing from one scan is often in the neighbours (the best charge SNR in the "
                     "window wins). 0 turns the rescue off, e.g. for native or targeted data.");
  defaults_.setMinInt("precursor_rescue_window", 0);
  defaults_.addTag("precursor_rescue_window", "advanced");
  defaults_.setValue("precursor_by_fragments", "true",
                     "The precursor of an MSn spectrum is chosen among the deconvolved target masses in its isolation window by the complementary "
                     "fragment pairs (b + y, c + z) of each mass beyond chance (the mean count at shifted masses): the spectrum is deconvolved once "
                     "with the bounds of all of them (the heaviest mass and the highest charge), and the mass of the highest charge SNR stays the "
                     "precursor unless another mass has at least two pairs beyond chance and two more than it; the spectrum is then deconvolved "
                     "with the chosen precursor's own bounds. With report_FDR, the decoy masses of the window then compete in the same way against "
                     "the chosen precursor, on the same fragment masses (with a single target mass in the window, those of its own deconvolution), "
                     "i.e., within the bounds of the target masses: a decoy mass heavier than all of them pairs only fragments up to the heaviest "
                     "one and can hardly win. A decoy mass that beats the precursor becomes the decoy precursor of the spectrum (FLASHTnT "
                     "-decoy_precursor); otherwise the decoy precursor is chosen among the window's decoy masses by the same rule (the one of the highest charge SNR unless the fragments single out another), so every MSn spectrum with a decoy mass in its window has one. With false, the precursor is the target mass and the decoy precursor the "
                     "decoy mass of the window of the highest charge SNR (unchecked: every MSn spectrum with a decoy mass in its window has one).");
  defaults_.setValidStrings("precursor_by_fragments", {"true", "false"});
  defaults_.addTag("precursor_by_fragments", "advanced");

  auto sd_defaults = SpectralDeconvolution().getDefaults();
  sd_defaults.remove("allowed_isotope_error");
  defaults_.insert("SD:", sd_defaults);

  Param mf_defaults = MassFeatureTrace().getDefaults();
  mf_defaults.setValue("min_cos", -1.0,
                       "Cosine similarity threshold between avg. and observed isotope pattern. When negative, MS1 cosine threshold for spectral "
                       "deconvolution (set by -SD:min_cos will be used ");
  mf_defaults.setValue(
    "mass_error_ppm", -1.0,
    "Specifies the mass error tolerance for feature tracing in ppm. When negative, the MS1 tolerance for deconvolution is used (e.g., 16 ppm is used when -SD:tol 16).");
  mf_defaults.addTag("min_cos", "advanced");
  mf_defaults.addTag("mass_error_ppm", "advanced");

  mf_defaults.remove("noise_threshold_int");
  mf_defaults.remove("reestimate_mt_sd");
  mf_defaults.remove("trace_termination_criterion");
  mf_defaults.remove("trace_termination_outliers");
  mf_defaults.remove("chrom_peak_snr");

  defaults_.insert("ft:", mf_defaults);
  defaults_.insert("iq:", TopDownIsobaricQuantification().getDefaults());

  defaultsToParam_();
}

void FLASHDeconvAlgorithm::updateMembers_()
{
  tols_ = param_.getValue("SD:tol");
  min_cos_ = param_.getValue("SD:min_cos");
  //precursor_MS1_window_ = param_.getValue("precursor_MS1_window");
  use_RNA_averagine_ = param_.getValue("use_RNA_averagine") != "false";
  report_decoy_ = param_.getValue("report_FDR") != "false";

  isolation_window_size_ = param_.getValue("isolation_window");
  scan_commands_file_ = param_.getValue("scan_commands").toString();
  precursor_rescue_window_ = param_.getValue("precursor_rescue_window");
  precursor_by_fragments_ = param_.getValue("precursor_by_fragments").toString() == "true";
  current_min_ms_level_ = max_ms_level_;
  current_max_ms_level_ = 0;
}

void FLASHDeconvAlgorithm::updateMSLevels_(MSExperiment& map)
{
  // read input dataset once to count spectra
  for (auto& it : map)
  {
    // if forced_ms_level > 0, force MS level of all spectra to 1.
    if (forced_ms_level_ > 0) { it.setMSLevel(forced_ms_level_); }

    if (it.empty()) { continue; }
    if (it.getMSLevel() > max_ms_level_) { continue; }

    uint ms_level = it.getMSLevel();
    current_max_ms_level_ = current_max_ms_level_ < ms_level ? ms_level : current_max_ms_level_;
    current_min_ms_level_ = current_min_ms_level_ > ms_level ? ms_level : current_min_ms_level_;
  }
  // Max MS Level is adjusted according to the input dataset
  current_max_ms_level_ = current_max_ms_level_ > max_ms_level_ ? max_ms_level_ : current_max_ms_level_;
}

void FLASHDeconvAlgorithm::filterLowPeaks_(MSExperiment& map)
{
  OPENMS_LOG_INFO << "Filtering low peaks in spectra ... ";
  ThresholdMower threshold_mower_filter;                         // threshold
  Param t_filter_param = threshold_mower_filter.getParameters(); //"threshold", .00001
  t_filter_param.setValue("threshold", 1e-6);
  threshold_mower_filter.setParameters(t_filter_param);

  // the threshold mower (filterPeakMap, spectrum by spectrum) and the trimming below in one parallel pass over the spectra, each thread
  // with its own copy of the mower; dynamic, as the spectra differ much in size
#pragma omp parallel default(none), shared(map, threshold_mower_filter)
  {
    ThresholdMower mower(threshold_mower_filter);
    std::vector<float> intensities; // scratch of the fast path
    PeakPickerHiRes picker;          // default parameters (no SNR filter)
#pragma omp for schedule(dynamic, 4)
    for (int i = 0; i < (int) map.size(); i++)
    {
      auto& it = map[i];
      // MSn profile spectra are centroided first (OpenMS PeakPickerHiRes, defaults). Deconvolved as profile, every profile point was a
      // peak: on 25 native top-down MS2 spectra (22,000-3,900,000 profile points) FLASHDeconv explained 605 fragments of the known
      // proteins (81 decoy-sequence matches) against 718 (78) after centroiding, as many as IsoDec, which centroids first, at the same
      // depth. MS1 profile spectra stay as they are (tolerance estimation, feature tracing and -quant use their profile points); the
      // raw MSn spectra of FLASHTnT's raw-spectrum scoring are copied before this step.
      if (it.getMSLevel() > 1 && it.getType(true) == SpectrumSettings::SpectrumType::PROFILE)
      {
        MSSpectrum picked;
        picker.pick(it, picked);
        it = std::move(picked);
      }
      mower.filterSpectrum(it);
      if (it.empty()) continue;
      Size count = it.getType(false) == SpectrumSettings::SpectrumType::CENTROID ? max_peak_count_for_centroid_ : max_peak_count_for_profile_;
      // The same result without the two sorts below: with strictly increasing m/z and no data arrays, the sort by position restores
      // the m/z order of the peaks the sort by intensity left, so the result is the peaks above the threshold in their order. The
      // threshold is the intensity of rank count (0-based) in descending order (none with at most count peaks) and the highest
      // intensity / 1000 (the mower left no NaN intensity).
      if (it.getFloatDataArrays().empty() && it.getStringDataArrays().empty() && it.getIntegerDataArrays().empty())
      {
        bool increasing = true;
        float max_intensity = it[0].getIntensity();
        for (Size k = 1; k < it.size(); k++)
        {
          if (! (it[k - 1].getMZ() < it[k].getMZ()))
          {
            increasing = false;
            break;
          }
          max_intensity = std::max(max_intensity, it[k].getIntensity());
        }
        if (increasing)
        {
          double threshold = 0;
          if (it.size() > count)
          {
            intensities.clear();
            for (const auto& p : it)
              intensities.push_back(p.getIntensity());
            std::nth_element(intensities.begin(), intensities.begin() + (SignedSize)count, intensities.end(), std::greater<float>());
            threshold = intensities[count];
          }
          threshold = std::max(threshold, (double)max_intensity / 1000);
          it.erase(std::remove_if(it.begin(), it.end(), [threshold](const Peak1D& p) { return p.getIntensity() <= threshold; }), it.end());
          continue;
        }
      }
      it.sortByIntensity(true);
      double threshold = it.size() <= count ? 0 : it[count].getIntensity();
      threshold = std::max(threshold, (double)it.begin()->getIntensity() / 1000);
      // pop back the low intensity peaks using threshold

      while (! it.empty() && it.back().getIntensity() <= threshold)
      {
        it.pop_back();
      }
      it.sortByPosition();
    }
  }
  OPENMS_LOG_INFO << "Done" << std::endl;
}

int FLASHDeconvAlgorithm::getScanNumber(const MSExperiment& map, Size index)
{
  auto native_id_str = map[index].getNativeID();
  std::vector<std::string> native_ids;
  StringUtils::split(native_id_str, ",", native_ids);
  std::string type_accession = "MS:1000768";

  if (!map.getSourceFiles().empty())
  {
    type_accession = map.getSourceFiles()[0].getNativeIDTypeAccession();
    if (type_accession.empty()) type_accession = "MS:1000768";
  }
  int scan_number = SpectrumLookup::extractScanNumber(native_ids.back(), type_accession);
  if (scan_number < 0) { scan_number = (int)index + 1; }

  return scan_number;
}

std::vector<double> FLASHDeconvAlgorithm::getTolerances() const
{
  return tols_;
}

struct FLASHDeconvAlgorithm::SpectrumJob_
{
  /// the spectrum: its index in the map and its scan number
  Size index = 0;
  int scan_number = 0;
  /// for an MSn spectrum without a precursor peak group: the scan number of the preceding spectrum of the lower MS level, if any
  bool has_fallback_scan = false;
  int fallback_scan = 0;

  /// round 1 ran (several target masses in the isolation window), and chose another precursor than the registered one (this one)
  bool several_window_masses = false;
  bool precursor_by_fragments = false;
  PeakGroup new_precursor_pg;
  Precursor new_precursor;
  /// round 2 ran (decoy masses in the window): the decoy precursor and whether a decoy mass beat the precursor
  bool window_decoys = false;
  bool decoy_by_fragments = false;
  double decoy_mass = 0, decoy_qscore = 0, decoy_qscore1d = -1;
  bool decoy_noise = false;

  /// the deconvolved spectrum, when it was deconvolved by a spectrum worker
  std::optional<DeconvolvedSpectrum> result;
};

void FLASHDeconvAlgorithm::deconvolveSpectrum_(const MSExperiment& map, SpectrumJob_& job, SpectralDeconvolution& sd,
                                               SpectralDeconvolution& noise_decoy_sd, SpectralDeconvolution& signal_decoy_sd) const
{
  const auto& spec = map[job.index];
  const uint ms_level = spec.getMSLevel();
  const int scan_number = job.scan_number;
  const std::string& native_id = spec.getNativeID();

  PeakGroup precursor_pg;
  const auto precursor_pg_it = native_id_precursor_peak_group_map_.find(native_id);
  if (precursor_pg_it != native_id_precursor_peak_group_map_.end()) precursor_pg = precursor_pg_it->second;

  // The precursor by the fragments, in two rounds on the same fragment masses. Round 1 (several target masses in the isolation
  // window): the MSn spectrum is deconvolved once within the bounds of all of them (the heaviest mass and the highest charge), the
  // target masses compete, and one whose complementary fragment pairs beat chance by the margin over the registered precursor's
  // replaces it. Round 2 (decoy masses in the window, i.e., report_FDR): the decoy masses of the window, as generated, compete against
  // the precursor of round 1 by the same score and margin on the fragment masses of round 1, or, with a single target mass in the
  // window, on those of its final deconvolution below (the same bounds, but for a rare window charge above that of the isolated ion);
  // decoy masses get no wider bounds. One that wins becomes the decoy precursor of the spectrum, and without a winner the spectrum has
  // none. Round 2 changes nothing of the target output. The spectrum is then deconvolved with the chosen precursor's own bounds (so
  // that the fragments reported stay within its mass and charge bounds). The precursor chosen by round 1 and the decoy precursor of
  // round 2 are registered for the spectrum by commitSpectrum_.
  auto sorted_masses = [](const DeconvolvedSpectrum& dspec) {
    std::vector<double> masses;
    masses.reserve(dspec.size());
    for (const auto& pg : dspec) masses.push_back(pg.getMonoMass());
    std::sort(masses.begin(), masses.end());
    return masses;
  };
  const std::vector<std::pair<double, bool>>* window_decoys = nullptr; // the challengers of round 2, if any
  std::vector<double> fragment_masses;                                   // the fragment masses both rounds score on
  bool round1 = false;                                                   // whether round 1 ran (several target masses)
  bool reuse_round1 = false; // round 1's deconvolution is the final one (same precursor, same mass and charge bounds)
  // not for a spectrum whose precursor FLASHIda commanded (-scan_commands): the commanded precursor stays (ADR-0046)
  if (precursor_by_fragments_ && ms_level > 1 && ! precursor_pg.empty() && located_by_scan_.find(scan_number) == located_by_scan_.end())
  {
    const auto cit = native_id_window_candidates_.find(native_id);
    const auto dit = native_id_window_decoys_.find(native_id);
    if (dit != native_id_window_decoys_.end() && ! dit->second.empty()) window_decoys = &dit->second;
    if (cit != native_id_window_candidates_.end() && cit->second.size() > 1)
    {
      round1 = true;
      job.several_window_masses = true;
      PeakGroup union_pg = precursor_pg;
      union_pg.setUseWindowBounds(true);
      sd.performSpectrumDeconvolution(spec, scan_number, union_pg);
      fragment_masses = sorted_masses(sd.getDeconvolvedSpectrum());

      std::vector<double> cand_masses;
      cand_masses.reserve(cit->second.size());
      for (const auto& cand : cit->second) cand_masses.push_back(cand.first.getMonoMass());
      const int best = choosePrecursorByFragments_(fragment_masses, cand_masses, precursor_pg.getMonoMass(), tols_[ms_level - 1]);
      if (best >= 0 && std::abs(cand_masses[best] - precursor_pg.getMonoMass()) > 1e-6)
      {
        job.precursor_by_fragments = true;
        precursor_pg = cit->second[best].first;
        job.new_precursor_pg = precursor_pg;
        job.new_precursor = cit->second[best].second;
      }
      else if ((int)sd.getParameters().getValue("precursor_charge") == 0 && (double)sd.getParameters().getValue("precursor_mz") == 0
               && (double)sd.getParameters().getValue("precursor_mass") == 0)
      {
        // The precursor stays. The deconvolution below differs from round 1 only by the window bounds; when they widen neither the mass
        // nor the charge bound (the heaviest and the most charged mass of the window are the precursor's), both are the same function of
        // the same inputs, and round 1's result is kept with the precursor as the deconvolution below would register it. Not with a
        // given precursor charge or m/z: the deconvolution then registers a precursor of its own instead of the one passed.
        auto& round1_spectrum = sd.getDeconvolvedSpectrum();
        const int union_charge = round1_spectrum.getCurrentMaxAbsCharge(std::numeric_limits<int>::max());
        const double union_mass = round1_spectrum.getCurrentMaxMass(std::numeric_limits<double>::max());
        round1_spectrum.setPrecursorPeakGroup(precursor_pg);
        reuse_round1 = round1_spectrum.getCurrentMaxAbsCharge(std::numeric_limits<int>::max()) == union_charge
                       && round1_spectrum.getCurrentMaxMass(std::numeric_limits<double>::max()) == union_mass;
      }
    }
  }
  // now do it
  if (! reuse_round1) sd.performSpectrumDeconvolution(spec, scan_number, precursor_pg);

  auto& deconvolved_spectrum = sd.getDeconvolvedSpectrum();

  // round 2, before the precursor information is added and the decoy masses of this spectrum are appended below
  if (window_decoys != nullptr)
  {
    job.window_decoys = true;
    if (! round1) fragment_masses = sorted_masses(deconvolved_spectrum);
    // the chosen precursor first (the default), then the decoy masses
    std::vector<double> cand_masses {precursor_pg.getMonoMass()};
    cand_masses.reserve(window_decoys->size() + 1);
    for (const auto& decoy : *window_decoys) cand_masses.push_back(decoy.first);
    int best = choosePrecursorByFragments_(fragment_masses, cand_masses, precursor_pg.getMonoMass(), tols_[ms_level - 1]);
    if (best > 0) { job.decoy_by_fragments = true; }
    else
    {
      // No decoy beats the precursor. The spectrum still gets a decoy precursor, chosen among the window's decoy masses by the rule
      // that chose the precursor among its target masses (the decoy of the highest charge SNR unless the fragments single out
      // another), so that every spectrum with a precursor from MS1 also has a decoy spectrum downstream (FLASHTnT -decoy_precursor):
      // with only the decoys that win here (1-3 of ~4100 spectra per E. coli lysate run), FLASHTnT saw almost no wrong precursor.
      std::vector<double> decoy_masses(cand_masses.begin() + 1, cand_masses.end());
      // (choosePrecursorByFragments_ returns -1 when the default, the first decoy, stays)
      best = 1 + std::max(0, choosePrecursorByFragments_(fragment_masses, decoy_masses, decoy_masses[0], tols_[ms_level - 1]));
    }
    job.decoy_mass = (*window_decoys)[best - 1].first;
    job.decoy_noise = (*window_decoys)[best - 1].second;
    // its precursor score as scored when the window was scanned
    job.decoy_qscore = 0;
    const auto qit = native_id_window_decoy_qscores_.find(native_id);
    if (qit != native_id_window_decoy_qscores_.end())
    {
      const auto q = qit->second.find(job.decoy_mass);
      if (q != qit->second.end()) job.decoy_qscore = q->second;
    }
    const auto q1it = native_id_window_decoy_qscores1d_.find(native_id);
    if (q1it != native_id_window_decoy_qscores1d_.end())
    {
      const auto q = q1it->second.find(job.decoy_mass);
      if (q != q1it->second.end()) job.decoy_qscore1d = q->second;
    }
  }

  // add precursor information
  if (precursor_pg_it != native_id_precursor_peak_group_map_.end())
  {
    if (job.precursor_by_fragments) { deconvolved_spectrum.setDeconvolvedPrecursor(job.new_precursor); }
    else
    {
      const auto pit = native_id_precursor_peak_map_.find(native_id);
      deconvolved_spectrum.setDeconvolvedPrecursor(pit != native_id_precursor_peak_map_.end() ? pit->second : Precursor());
    }
    // the decoy precursor: of round 2, else as registered with the precursor (without precursor_by_fragments)
    bool has_decoy = job.window_decoys;
    double decoy_mass = job.decoy_mass, decoy_qscore = job.decoy_qscore, decoy_qscore1d = job.decoy_qscore1d;
    bool decoy_noise = job.decoy_noise;
    if (! has_decoy)
    {
      const auto decoy_it = native_id_decoy_precursor_mass_map_.find(native_id);
      if (decoy_it != native_id_decoy_precursor_mass_map_.end())
      {
        has_decoy = true;
        decoy_mass = decoy_it->second;
        const auto qit = native_id_decoy_precursor_qscore_map_.find(native_id);
        decoy_qscore = qit != native_id_decoy_precursor_qscore_map_.end() ? qit->second : 0;
        const auto q1it = native_id_decoy_precursor_qscore1d_map_.find(native_id);
        decoy_qscore1d = q1it != native_id_decoy_precursor_qscore1d_map_.end() ? q1it->second : -1;
        const auto nit = native_id_decoy_precursor_noise_map_.find(native_id);
        decoy_noise = nit != native_id_decoy_precursor_noise_map_.end() && nit->second;
      }
    }
    if (has_decoy)
    {
      deconvolved_spectrum.setDecoyPrecursorMass(decoy_mass);
      deconvolved_spectrum.setDecoyPrecursorQscore(decoy_qscore);
      deconvolved_spectrum.setDecoyPrecursorQscore1D(decoy_qscore1d);
      // a noise decoy counts with the noise decoy weight, known after the q-values below
      deconvolved_spectrum.setDecoyPrecursorWeight(decoy_noise ? -1.0 : 1.0);
    }
  }

  const auto mit = native_id_precursor_mass_intensity_map_.find(native_id);
  if (mit != native_id_precursor_mass_intensity_map_.end()) { deconvolved_spectrum.setPrecursorMassIntensityMap(mit->second); }
  // add precursor scan number information when a precursor peak group is not found
  if (ms_level > 1 && precursor_pg.empty() && job.has_fallback_scan) { deconvolved_spectrum.setPrecursorScanNumber(job.fallback_scan); }

  if (report_decoy_ && ! deconvolved_spectrum.empty())
  {
    // the two decoy runs in parallel, but one after the other in a spectrum worker
    bool decoy_sections = true;
#ifdef _OPENMP
    decoy_sections = ! omp_in_parallel();
#endif
    if (decoy_sections)
    {
#pragma omp parallel sections default(none) shared(spec, scan_number, precursor_pg, noise_decoy_sd, signal_decoy_sd)
      {
#pragma omp section
        noise_decoy_sd.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
#pragma omp section
        signal_decoy_sd.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
      }
    }
    else
    {
      noise_decoy_sd.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
      signal_decoy_sd.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
    }
    deconvolved_spectrum.sortByQscore();

    deconvolved_spectrum.reserve(deconvolved_spectrum.size() + signal_decoy_sd.getDeconvolvedSpectrum().size()
                                 + noise_decoy_sd.getDeconvolvedSpectrum().size());

    // the decoy spectra are rebuilt for the next spectrum, so their peak groups are moved
    for (auto& pg : signal_decoy_sd.getDeconvolvedSpectrum())
      deconvolved_spectrum.push_back(std::move(pg));

    for (auto& pg : noise_decoy_sd.getDeconvolvedSpectrum())
      deconvolved_spectrum.push_back(std::move(pg));
    noise_decoy_sd.releasePools();
    signal_decoy_sd.releasePools();

    deconvolved_spectrum.sort();
  }
}

void FLASHDeconvAlgorithm::commitSpectrum_(const MSExperiment& map, const SpectrumJob_& job)
{
  const std::string& native_id = map[job.index].getNativeID();
  if (job.several_window_masses) spectra_with_several_window_masses_++;
  if (job.precursor_by_fragments)
  {
    precursors_by_fragments_++;
    native_id_precursor_peak_group_map_[native_id] = job.new_precursor_pg;
    native_id_precursor_peak_map_[native_id] = job.new_precursor;
  }
  else if (native_id_precursor_peak_group_map_.find(native_id) != native_id_precursor_peak_group_map_.end())
  {
    native_id_precursor_peak_map_[native_id]; // a precursor peak group without a precursor peak (FLASHIda log) gets an empty one, as before
  }
  if (job.window_decoys)
  {
    spectra_with_window_decoys_++;
    if (job.decoy_by_fragments) decoy_precursors_by_fragments_++;
    decoy_precursors_++;
    native_id_decoy_precursor_mass_map_[native_id] = job.decoy_mass;
    native_id_decoy_precursor_noise_map_[native_id] = job.decoy_noise;
    native_id_decoy_precursor_qscore_map_[native_id] = job.decoy_qscore;
    native_id_decoy_precursor_qscore1d_map_[native_id] = job.decoy_qscore1d;
  }
}

void FLASHDeconvAlgorithm::runSpectralDeconvolution_(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
{
  startProgress(0, (SignedSize)map.size(), "running FLASHDeconv");
  std::vector<int> scan_numbers(map.size());
  for (Size index = 0; index < map.size(); index++)
  {
    int scan_number = getScanNumber(map, index);
    scan_numbers[index] = scan_number;
  }

  // Spectrum workers: the spectra of an MS level are deconvolved in parallel, whole spectra per thread, each thread with its own
  // SpectralDeconvolution instances (target, and with report_FDR the noise and signal decoys) sharing one immutable averagine; the
  // per-peak-group loops of SpectralDeconvolution then run serially inside a worker. The spectra only read the precursor maps (filled
  // per MS level before, by findPrecursorPeakGroupsForMSnSpectra_), and what they change (the precursor chosen by the fragments, the decoy
  // precursor, the counters) is committed after the join in the input order, as are the deconvolved spectra: the output is that of the
  // serial run. Not with a FLASHIda log, a given precursor charge or m/z, inside a parallel
  // region, with a single thread, or for an MS level with repeated native IDs (a spectrum would read what another one changes): there the
  // spectra are deconvolved one after the other.
  struct SpectrumWorker
  {
    SpectralDeconvolution sd, noise_decoy_sd, signal_decoy_sd;
  };
  std::vector<SpectrumWorker> workers;
#ifdef _OPENMP
  const Param sd_param = sd_.getParameters();
  if (scan_commands_file_.empty() && ! omp_in_parallel() && (int)sd_param.getValue("precursor_charge") == 0
      && (double)sd_param.getValue("precursor_mz") == 0 && (double)sd_param.getValue("precursor_mass") == 0)
  {
    std::vector<Size> level_counts(current_max_ms_level_ + 1, 0);
    for (const auto& spec : map)
    {
      if (! spec.empty() && spec.getMSLevel() >= 1 && spec.getMSLevel() <= current_max_ms_level_) level_counts[spec.getMSLevel()]++;
    }
    const int worker_count = (int)std::min((Size)omp_get_max_threads(), *std::max_element(level_counts.begin(), level_counts.end()));
    if (worker_count > 1)
    {
      // constructed in place once: the decoy instances keep a pointer to the deconvolved spectrum of their worker's target instance
      workers = std::vector<SpectrumWorker>(worker_count);
      const auto avg = std::make_shared<const FLASHHelperClasses::PrecalculatedAveragine>(sd_.getAveragine());
      std::shared_ptr<const FLASHHelperClasses::PrecalculatedAveragine> noise_decoy_avg;
      if (report_decoy_) noise_decoy_avg = std::make_shared<const FLASHHelperClasses::PrecalculatedAveragine>(sd_noise_decoy_.getAveragine());
      for (auto& worker : workers)
      {
        worker.sd.setParameters(sd_param);
        worker.sd.setSharedAveragine_(avg);
        if (! report_decoy_) continue;
        worker.noise_decoy_sd.setParameters(sd_noise_decoy_.getParameters());
        worker.noise_decoy_sd.setTargetDecoyType(PeakGroup::TargetDecoyType::noise_decoy, worker.sd.getDeconvolvedSpectrum());
        worker.noise_decoy_sd.setSharedAveragine_(noise_decoy_avg);
        worker.signal_decoy_sd.setParameters(sd_signal_decoy_.getParameters());
        worker.signal_decoy_sd.setTargetDecoyType(PeakGroup::TargetDecoyType::signal_decoy, worker.sd.getDeconvolvedSpectrum());
        worker.signal_decoy_sd.setSharedAveragine_(avg);
      }
    }
  }
#endif

  for (uint ms_level = 1; ms_level <= current_max_ms_level_; ms_level++)
  {
    if (ms_level > 1)
    {
      // here, register precursor peak groups to the ms2 spectra.
      findPrecursorPeakGroupsForMSnSpectra_(map, deconvolved_spectra, ms_level);
    }

    // the precursor scan of an MSn spectrum without a precursor peak group: the preceding spectrum of the lower MS level
    auto needs_fallback_scan = [&](const MSSpectrum& spec) {
      if (ms_level <= 1) return false;
      const auto it = native_id_precursor_peak_group_map_.find(spec.getNativeID());
      return it == native_id_precursor_peak_group_map_.end() || it->second.empty();
    };

    if (! workers.empty())
    {
      std::vector<SpectrumJob_> jobs;
      std::set<std::string> native_ids;
      bool unique_native_ids = true;
      Size progress_count = 0;
      Size lower_index = map.size(); // the last spectrum of the lower MS level so far
      for (Size index = 0; index < map.size(); index++)
      {
        const auto& spec = map[index];
        if (spec.getMSLevel() == ms_level - 1) lower_index = index;
        if (spec.getMSLevel() != ms_level) { continue; }
        progress_count++;
        if (spec.empty()) { continue; }
        if (! native_ids.insert(spec.getNativeID()).second)
        {
          unique_native_ids = false;
          break;
        }
        SpectrumJob_ job;
        job.index = index;
        job.scan_number = scan_numbers[index];
        if (lower_index < map.size() && needs_fallback_scan(spec))
        {
          job.has_fallback_scan = true;
          job.fallback_scan = scan_numbers[lower_index];
        }
        jobs.push_back(std::move(job));
      }

      if (unique_native_ids && jobs.size() > 1)
      {
        std::vector<std::exception_ptr> errors(jobs.size());
        const int team_size = (int)std::min(workers.size(), jobs.size());
        // the most expensive spectra (peak count as the cost proxy) start first, which shortens the tail of the wave; the results stay in
        // their job slots, so the output order is unchanged
        std::vector<SignedSize> job_order(jobs.size());
        std::iota(job_order.begin(), job_order.end(), 0);
        std::stable_sort(job_order.begin(), job_order.end(), [&](SignedSize a, SignedSize b) { return map[jobs[a].index].size() > map[jobs[b].index].size(); });
        // The features of the minimum MS level depend only on its deconvolved spectra (all of deconvolved_spectra here), which the spectra
        // of this level neither read nor change: when this is the last level, one thread traces them (runFeatureFinding_ takes them) while
        // the others deconvolve. The tracing visits the spectra in the order the serial path does (in the scan order of all spectra) only
        // when their scan numbers are unique, so only then.
        bool trace_in_wave = ms_level == current_min_ms_level_ + 1 && ms_level == current_max_ms_level_;
        for (Size i = 1; trace_in_wave && i < deconvolved_spectra.size(); i++)
        {
          if (deconvolved_spectra[i].getScanNumber() == deconvolved_spectra[i - 1].getScanNumber()) trace_in_wave = false;
        }
#pragma omp parallel num_threads(team_size) default(none) shared(map, jobs, workers, errors, job_order, deconvolved_spectra, trace_in_wave)
        {
          int worker_index = 0;
#ifdef _OPENMP
          worker_index = omp_get_thread_num();
#endif
          auto& worker = workers[worker_index];
          if (trace_in_wave)
          {
#pragma omp single nowait
            {
              try
              {
                min_level_features_ = traceMinLevelFeatures_(deconvolved_spectra);
              }
              catch (...)
              {
                min_level_trace_error_ = std::current_exception();
              }
            }
          }
// the costliest spectra first (job_order) only with a monotonic schedule (nonmonotonic is the OpenMP >= 5.0 default, e.g. clang/libomp)
#if defined(_OPENMP) && _OPENMP >= 201511
  #pragma omp for schedule(monotonic: dynamic, 1)
#else
  #pragma omp for schedule(dynamic, 1)
#endif
          for (SignedSize position = 0; position < (SignedSize)jobs.size(); position++)
          {
            const SignedSize j = job_order[position];
            try
            {
              deconvolveSpectrum_(map, jobs[j], worker.sd, worker.noise_decoy_sd, worker.signal_decoy_sd);
              jobs[j].result.emplace(std::move(worker.sd.getDeconvolvedSpectrum()));
            }
            catch (...)
            {
              errors[j] = std::current_exception();
            }
          }
        }
        for (const auto& error : errors)
        {
          if (! error) continue;
          endProgress();
          std::rethrow_exception(error);
        }
        // commit in the input order
        deconvolved_spectra.reserve(deconvolved_spectra.size() + jobs.size());
        for (auto& job : jobs)
        {
          commitSpectrum_(map, job);
          deconvolved_spectra.push_back(std::move(*job.result));
        }
        for (Size i = 0; i < progress_count; i++)
          nextProgress();
        std::sort(deconvolved_spectra.begin(), deconvolved_spectra.end());
        continue;
      }
    }

    // run Spectral deconvolution, one spectrum after the other
    Size lower_index = map.size(); // the last spectrum of the lower MS level so far
    for (Size index = 0; index < map.size(); index++)
    {
      const auto& spec = map[index];
      if (spec.getMSLevel() == ms_level - 1) lower_index = index;
      if (ms_level != spec.getMSLevel()) { continue; }

      nextProgress();
      if (spec.empty()) { continue; }

      SpectrumJob_ job;
      job.index = index;
      job.scan_number = getScanNumber(map, index);
      if (lower_index < map.size() && needs_fallback_scan(spec))
      {
        job.has_fallback_scan = true;
        job.fallback_scan = getScanNumber(map, lower_index);
      }
      deconvolveSpectrum_(map, job, sd_, sd_noise_decoy_, sd_signal_decoy_);
      commitSpectrum_(map, job);
      deconvolved_spectra.push_back(std::move(sd_.getDeconvolvedSpectrum()));
    }
    std::sort(deconvolved_spectra.begin(), deconvolved_spectra.end());
  }
  endProgress();
  // each part only if its case occurred in this run (decoy masses in a window only with report_FDR)
  if (spectra_with_several_window_masses_ > 0 || spectra_with_window_decoys_ > 0)
  {
    std::string by_fragments = "Precursors by the fragments", sep = ": ";
    if (spectra_with_several_window_masses_ > 0)
    {
      by_fragments += sep + "masses chosen among the target masses of the isolation window in " + StringUtils::toStr(precursors_by_fragments_) + " of "
                      + StringUtils::toStr(spectra_with_several_window_masses_) + " MSn spectra with several";
      sep = "; ";
    }
    if (spectra_with_window_decoys_ > 0)
    {
      by_fragments += sep + "decoy precursors in " + StringUtils::toStr(decoy_precursors_) + " (a decoy mass of the window beat the chosen precursor in " + StringUtils::toStr(decoy_precursors_by_fragments_) + ")"
                      + " of " + StringUtils::toStr(spectra_with_window_decoys_) + " MSn spectra with a precursor and decoy masses in the window";
    }
    OPENMS_LOG_INFO << by_fragments << std::endl;
  }
}


const FLASHHelperClasses::PrecalculatedAveragine& FLASHDeconvAlgorithm::getAveragine()
{
  return sd_.getAveragine();
}

const FLASHHelperClasses::PrecalculatedAveragine& FLASHDeconvAlgorithm::getDecoyAveragine()
{
  if (!report_decoy_) return getAveragine();
  return sd_noise_decoy_.getAveragine();
}

namespace
{
  /// The half-width (ppm) holding @p coverage of the true peaks' errors among @p abs_errors (|ppm| of the peaks of good masses
  /// recruited within @p max_error): the errors histogrammed in 0.1 ppm bins, the background of random peaks taken as the mean
  /// count over [0.5, 0.75) x max_error (beyond the true peaks' tails, below the rise at the window's edge), and the excess over it
  /// summed from 0. -1 if there is no excess or the half-width reaches 0.5 x max_error (the window was too narrow).
  /// @p centre_excess: the mean excess per bin within the half-width relative to the background (how far the true peaks stand out of
  /// the random ones; near 0 when the true errors are spread flat over the whole window and the half-width is noise)
  double trueErrorHalfWidth(const std::vector<double>& abs_errors, const double max_error, const double coverage, double& centre_excess)
  {
    centre_excess = 0;
    const double bin = 0.1;
    const Size bin_count = (Size)std::round(max_error / bin);
    std::vector<double> counts(bin_count, 0);
    for (const double e : abs_errors)
    {
      if (e < max_error) counts[std::min(bin_count - 1, (Size)(e / bin))]++;
    }
    const Size band_begin = bin_count / 2, band_end = bin_count * 3 / 4;
    double background = 0;
    for (Size b = band_begin; b < band_end; b++) background += counts[b];
    background /= (double)(band_end - band_begin);
    double total = 0;
    for (Size b = 0; b < band_begin; b++) total += counts[b] - background;
    if (total <= 0) return -1;
    double sum = 0;
    for (Size b = 0; b < band_begin; b++)
    {
      sum += counts[b] - background;
      if (sum < coverage * total) continue;
      const double centre = sum / (double)(b + 1); // the mean excess per bin within the half-width
      centre_excess = background > 0 ? centre / background : std::numeric_limits<double>::infinity();
      return (double)(b + 1) * bin;
    }
    return -1;
  }
} // namespace

void FLASHDeconvAlgorithm::determineTolerance_(const MSExperiment& map, const Param& sd_param, const FLASHHelperClasses::PrecalculatedAveragine& avg, const uint ms_level)
{
  // The tolerance that recruits 98% of the true peaks: a peak is taken when its error from its mass's isotope position is within
  // +-0.8 x tol (PeakGroup), so tol = W98 / 0.8, W98 the half-width holding 98% of the errors of the peaks of good masses (Qscore >= 0.9)
  // deconvolved with a wide tolerance, after the random peaks recruited with them are subtracted. The errors are relative to each
  // mass's own fitted mass, so a calibration offset does not show (nor does it affect the deconvolution; the median is logged).
  // The coverage, not the spread: the errors are heavy-tailed (the 98% half-width is 2.7 to 6.3 times the robust sigma), so a multiple
  // of sigma gave too narrow a tolerance on some data. Eclip_1480 MS2 (FLASHTnT one step, proteoforms at 1%): 5 ppm 1565, 7 ppm 1742,
  // 10 ppm 1716, 15 ppm 1692; this rule gives 8.0. Replaces a Gaussian fit to 1-ppm bins at 200 ppm, which failed on most MS1 data.
  OPENMS_LOG_INFO << "Determining tolerance for MS" << ms_level << " ... " << std::flush;
  const auto start_time = std::chrono::steady_clock::now();
  std::vector<Size> spectra; // the non-empty spectra of the level, in file order
  for (Size i = 0; i < map.size(); i++)
  {
    if ((uint)map[i].getMSLevel() == ms_level && ! map[i].empty()) spectra.push_back(i);
  }
  if (spectra.empty())
  {
    OPENMS_LOG_INFO << "failed. Cannot be determined - no MS" << ms_level << " spectrum. Set to 10ppm (default tolerance)." << std::endl;
    tols_[ms_level - 1] = 10;
    return;
  }
  const bool centroid = map[spectra[0]].getType(false) != SpectrumSettings::SpectrumType::PROFILE;
  const auto shared_avg = std::make_shared<const FLASHHelperClasses::PrecalculatedAveragine>(avg);
  const int min_charge = sd_param.getValue("min_charge");

  const double coverage = .98, recruit = .8, min_centre_excess = .25;
  // the errors of the peaks of good masses (or, for profile spectra, per charge and isotope slot, the intensity-weighted mean) of the
  // spectra at @p positions, deconvolved with tolerance @p sample_tol, in parallel; one vector per position, so the result does not
  // depend on the number of threads. In a window wider than 30 ppm, only masses whose isotope spacing (ppm: the isotope distance / mass)
  // is at least twice the recruitment window: in a heavier mass a peak of the next isotope falls in the window, and its errors spread to
  // half the spacing instead of to the window's edge (no background band; at 300 ppm, MS1 of an ETD file with 5 ppm added to every peak
  // gave 64 ppm). Not at 30 ppm: there the heavier masses' wider errors are those FD meets at its working tolerances (all masses: 1480
  // MS1 6.6 ppm; masses below 20.8 kDa only: 5.1 ppm, and MS1 5 ppm lost proteoforms against 10 ppm in FLASHTnT, 1688 vs 1716).
  const auto collect = [&](const std::vector<Size>& positions, const double sample_tol) {
    // in isotope distances (x getIsotopeDaDistance())
    const double max_mass = sample_tol > 30.0 ? 1e6 / (2.0 * recruit * sample_tol) : std::numeric_limits<double>::max();
    Param p = sd_param;
    DoubleList tols(tols_.size(), sample_tol);
    p.setValue("tol", tols);
    p.setValue("min_charge", min_charge < 0 ? -1 : 1); // charge 1 and small masses give many peaks to measure
    p.setValue("min_mass", 50.0);
    std::vector<std::vector<double>> errors(positions.size());
    std::vector<std::exception_ptr> failures(positions.size());
#pragma omp parallel
    {
      SpectralDeconvolution sd;
      sd.setParameters(p);
      sd.setSharedAveragine_(shared_avg);
#pragma omp for schedule(dynamic, 1)
      for (SignedSize k = 0; k < (SignedSize)positions.size(); k++)
      {
        try
        {
          sd.performSpectrumDeconvolution(map[spectra[positions[k]]], 0, PeakGroup());
          for (const auto& pg : sd.getDeconvolvedSpectrum())
          {
            if (pg.getQscore() < .9 || pg.getMonoMass() > max_mass * pg.getIsotopeDaDistance()) continue;
            std::map<std::pair<int, int>, std::pair<double, double>> per_slot; // (charge, isotope) -> (sum of intensity x error, sum of intensity)
            for (const auto& peak : pg)
            {
              const double e = ((pg.getMonoMass() + peak.isotopeIndex * pg.getIsotopeDaDistance()) / peak.abs_charge
                                + FLASHHelperClasses::getChargeMass(peak.is_positive) - peak.mz) / peak.mz * 1e6;
              if (centroid) { errors[k].push_back(e); continue; }
              auto& slot = per_slot[{peak.abs_charge, peak.isotopeIndex}];
              slot.first += peak.intensity * e;
              slot.second += peak.intensity;
            }
            for (const auto& [key, slot] : per_slot)
            {
              if (slot.second > 0) errors[k].push_back(slot.first / slot.second);
            }
          }
        }
        catch (...)
        {
          failures[k] = std::current_exception();
        }
      }
    }
    for (const auto& failure : failures)
    {
      if (failure) std::rethrow_exception(failure);
    }
    return errors;
  };

  // The sample tolerance (window) widens until the true errors end well inside it: at 30 ppm the background band starts at 12 ppm, so
  // true errors reaching it raise the background and cap the half-width near 11 ppm (tol ~14) however wide they are. A window is taken
  // when its half-width is at most a third of the recruitment window (the true peaks end before the band) and the true peaks stand out
  // (centre excess at least 0.25 x background; a flat spread gives a half-width of noise); else the widest window with a half-width whose
  // peaks stand out (it may underestimate); else 10 ppm. Simulated Gaussian and t(3) errors, random peaks 15% to 150% of the true ones:
  // tol needed up to ~40 ppm is found within 15%; 30 ppm alone gave 14 for every tol needed above 14.
  // Per window, spectra are sampled in rounds, each doubling the sample, until 20,000 errors or every 4th spectrum.
  const std::vector<std::pair<Size, Size>> rounds {{64, 0}, {64, 32}, {32, 16}, {16, 8}, {8, 4}};
  const std::vector<double> sample_tols {30.0, 100.0, 300.0};
  const std::string unit = centroid ? "peak errors" : "slot errors (profile)";
  std::vector<double> sampled;   // the errors of the window taken
  std::string trace;                  // the windows tried, for the log
  Size sampled_spectra = 0, max_error_count = 0;
  double half_width = -1, sample_tol = 0, centre_excess = 0;
  bool within_window = false; // the window taken holds the true errors well inside (not only the widest with a half-width)
  for (const double tol : sample_tols)
  {
    std::vector<double> errors;
    Size spectrum_count = 0;
    for (const auto& [stride, offset] : rounds)
    {
      std::vector<Size> positions;
      for (Size i = offset; i < spectra.size(); i += stride) positions.push_back(i);
      spectrum_count += positions.size();
      for (const auto& e : collect(positions, tol)) errors.insert(errors.end(), e.begin(), e.end());
      if (errors.size() >= 20000) break;
    }
    max_error_count = std::max(max_error_count, errors.size());
    trace += std::string(trace.empty() ? "" : ", ") + StringUtils::toStr((int)tol) + " ppm: " + StringUtils::toStr(errors.size()) + " errors";
    if (errors.size() < 1000) continue;
    std::vector<double> abs_errors(errors.size());
    std::transform(errors.begin(), errors.end(), abs_errors.begin(), [](double e) { return std::abs(e); });
    double excess = 0;
    const double hw = trueErrorHalfWidth(abs_errors, recruit * tol, coverage, excess);
    trace += ", half-width " + (hw > 0 ? StringUtils::number(hw, 1) : std::string("none")) + ", centre " + StringUtils::number(excess, 1) + " x background";
    // no random peak in the band (an infinite centre excess) means the window is not measured beyond the true peaks: not taken
    if (hw <= 0 || excess < min_centre_excess || std::isinf(excess)) continue;
    half_width = hw;
    sample_tol = tol;
    centre_excess = excess;
    sampled = std::move(errors);
    sampled_spectra = spectrum_count;
    within_window = hw <= recruit * tol / 3.0;
    if (within_window) break;
  }
  if (half_width < 0)
  {
    if (max_error_count < 1000)
    {
      OPENMS_LOG_INFO << "failed. Only " << max_error_count << " " << unit << " from Qscore >= 0.9 masses (1000 needed; " << trace
                      << "). Set to 10ppm (default tolerance)." << std::endl;
    }
    else
    {
      OPENMS_LOG_INFO << "failed. The " << unit << " show no true-peak distribution standing out of the random peaks (" << trace
                      << "). Set to 10ppm (default tolerance)." << std::endl;
    }
    tols_[ms_level - 1] = 10;
    return;
  }
  tols_[ms_level - 1] = std::round(10.0 * half_width / recruit) / 10.0;
  // the centre and the robust sigma, logged only
  std::vector<double> sorted = sampled;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  for (auto& e : sorted) e = std::abs(e - median);
  std::sort(sorted.begin(), sorted.end());
  const double sigma = 1.4826 * sorted[sorted.size() / 2];
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
  OPENMS_LOG_INFO << "done. Determined tolerance: " << tols_[ms_level - 1] << " ppm (" << (int)std::round(100 * coverage) << "% of " << sampled.size() << " " << unit
                  << " within +-" << StringUtils::number(half_width, 1) << " ppm = " << recruit << " x tol; " << sampled_spectra << " of " << spectra.size()
                  << " spectra deconvolved at " << sample_tol << " ppm" << (within_window ? "" : ", the widest window: may be underestimated")
                  << "; true peaks " << StringUtils::number(centre_excess, 1) << " x background at the centre; centre " << StringUtils::number(median, 2)
                  << " ppm, robust sigma " << StringUtils::number(sigma, 2) << " ppm; " << StringUtils::number(seconds, 1) << " s; windows " << trace << ")." << std::endl;
}

void FLASHDeconvAlgorithm::run(MSExperiment& map,
                               std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                               std::vector<FLASHHelperClasses::MassFeature>& deconvolved_features)
{
  min_level_features_.reset();
  min_level_trace_error_ = nullptr;
  MassFeatureTrace::resetFeatureIndex(); // the feature indices of this run start at 1

  // ADR-0046: join this data file to its acquisition BEFORE any deconvolution, so a scan_commands file that
  // does not belong to it aborts in seconds. Metadata only -- nothing below has touched the spectra yet.
  located_by_scan_.clear();
  if (!scan_commands_file_.empty())
  {
    std::vector<ScanCommandJoin::Scan> scans;
    scans.reserve(map.size());
    for (Size index = 0; index < map.size(); index++)
    {
      ScanCommandJoin::Scan s;
      s.scan_number = getScanNumber(map, index);
      s.ms_level = (int)map[index].getMSLevel();
      if (map[index].metaValueExists("scan description")) { s.scan_description = map[index].getMetaValue("scan description").toString(); }
      for (const auto& precursor : map[index].getPrecursors())
      {
        s.isolation_targets.push_back(precursor.getMZ());
        // MzMLHandler overwrites getMZ() with the selected-ion m/z when the two differ, and keeps the target here
        if (precursor.metaValueExists("isolation window target m/z")) { s.isolation_targets.push_back((double)precursor.getMetaValue("isolation window target m/z")); }
      }
      scans.push_back(s);
    }
    located_by_scan_ = ScanCommandJoin::join(scans, ScanCommandJoin::parse(scan_commands_file_));
    OPENMS_LOG_INFO << "FLASHIda scan commands used: " << scan_commands_file_ << " (" << located_by_scan_.size() << " commanded MS2 scans located)" << std::endl;
  }

  updateMSLevels_(map);
  filterLowPeaks_(map);

  sd_ = SpectralDeconvolution();
  Param sd_param = param_.copy("SD:", true);
  sd_param.setValue("allowed_isotope_error", param_.getValue("allowed_isotope_error"));
  OPENMS_LOG_INFO<< "Calculating Averagines ... " << std::flush;
  sd_.setParameters(sd_param);
  sd_.calculateAveragine(use_RNA_averagine_);
  OPENMS_LOG_INFO<< "Done" << std::endl;
  const auto& avg = sd_.getAveragine();

  // determine tolerance in case tolerance input is negative
  for (uint ms_level = 1; ms_level <= current_max_ms_level_; ms_level++)
  {
    if (tols_[ms_level - 1] > 0) continue;
    determineTolerance_(map, sd_param, avg, ms_level);
  }

  sd_param.setValue("tol", tols_);
  sd_.setParameters(sd_param);

  if (report_decoy_)
  {
    sd_noise_decoy_.setParameters(sd_param);
    sd_noise_decoy_.setTargetDecoyType(PeakGroup::TargetDecoyType::noise_decoy, sd_.getDeconvolvedSpectrum()); // noise
    sd_noise_decoy_.calculateAveragine(use_RNA_averagine_); // for noise, averagine needs to be calculated differently.

    sd_signal_decoy_.setParameters(sd_param);
    sd_signal_decoy_.setTargetDecoyType(PeakGroup::TargetDecoyType::signal_decoy, sd_.getDeconvolvedSpectrum()); // isotope
    sd_signal_decoy_.setAveragine(avg);
  }

  setLogType(CMD);
  deconvolved_spectra.reserve(map.size() * 4);

  // run spectral deconvolution here and get deconvolved spectra
  runSpectralDeconvolution_(map, deconvolved_spectra);

  // feature tracing here and update FeatureQScores
  runFeatureFinding_(deconvolved_spectra, deconvolved_features);

  noise_decoy_weights_ = Qvalue::updatePeakGroupQvalues(deconvolved_spectra);
  if (report_decoy_)
  {
    std::string weights;
    for (const auto& [level, weight] : noise_decoy_weights_) weights += " MS" + StringUtils::toStr(level) + " " + StringUtils::number(weight, 3);
    OPENMS_LOG_INFO << "Noise decoy weights:" << weights << std::endl;
  }
  for (auto& dspec : deconvolved_spectra)
  {
    // a decoy precursor is a mass of the level below: its weight is that level's
    const uint ms_level = dspec.getOriginalSpectrum().getMSLevel();
    if (dspec.getDecoyPrecursorMass() > 0 && dspec.getDecoyPrecursorWeight() < 0) dspec.setDecoyPrecursorWeight(getNoiseDecoyWeight(ms_level > 1 ? ms_level - 1 : 1));
  }

  TopDownIsobaricQuantification quantifier;
  Param quant_param = param_.copy("iq:", true);
  quantifier.setParameters(quant_param);
  // Isobaric quant run
  quantifier.quantify(map, deconvolved_spectra, deconvolved_features);
}

// ADR-0046 decision 3: FLASHDeconv's deconvolution of the commanded survey holds NOTHING in the isolation window,
// so the precursor is rebuilt from the command itself -- one peak at the anchor charge, scores as the engine
// logged them. The construction the ida.log path used, over doubles instead of floats, and with no charge range
// to claim. A window that holds something the command did not name takes the charge-SNR pick instead: a measured
// species carries an isotope envelope, and this one does not -- which also bounds the FRAGMENT deconvolution
// (SpectralDeconvolution.cpp:197-218), so a rebuilt precursor costs the MS2 as well as the mass.
PeakGroup FLASHDeconvAlgorithm::peakGroupFromCommand_(const ScanCommandJoin::Located& located) const
{
  const auto& row = located.row;
  FLASHHelperClasses::LogMzPeak precursor_log_mz_peak;
  precursor_log_mz_peak.abs_charge = std::abs(row.charge);
  precursor_log_mz_peak.is_positive = row.charge > 0;
  precursor_log_mz_peak.isotopeIndex = 0;
  precursor_log_mz_peak.mass = row.mono_mass;
  precursor_log_mz_peak.intensity = (float)row.peakgroup_intensity;

  PeakGroup precursor_pg(precursor_log_mz_peak.abs_charge, precursor_log_mz_peak.abs_charge, precursor_log_mz_peak.is_positive);
  precursor_pg.push_back(precursor_log_mz_peak);
  precursor_pg.setChargeIsotopeCosine(precursor_log_mz_peak.abs_charge, (float)row.charge_cos);
  precursor_pg.setChargeSNR(precursor_log_mz_peak.abs_charge, (float)row.charge_snr);
  precursor_pg.setIsotopeCosine((float)row.iso_cos);
  precursor_pg.setSNR((float)row.snr);
  precursor_pg.setChargeScore((float)row.charge_score);
  precursor_pg.setAvgPPMError((float)row.ppm_error);
  precursor_pg.setQscore(row.qscore);
  precursor_pg.setRepAbsCharge(precursor_log_mz_peak.abs_charge);
  precursor_pg.updateMonoMassAndIsotopeIntensities(tols_[0]);
  precursor_pg.setScanNumber(located.parent_scan_number);
  return precursor_pg;
}

int FLASHDeconvAlgorithm::choosePrecursorByFragments_(const std::vector<double>& fragment_masses, const std::vector<double>& cand_masses,
                                                      double default_mass, double tol_ppm)
{
  // complementary pairs summing to the mass (b + y) or the mass + H (c + z) within twice the tolerance. Two masses less than 10 Da apart
  // are no pair: a mass at about half the candidate (e.g., the unfragmented precursor deconvolved at half its mass, with its isotope and
  // hydrogen variants) would pair with itself
  auto pairs_at = [&](double mass) {
    int n = 0;
    for (double band : {0.0, Constants::PROTON_MASS_U})
    {
      const double target = mass + band, tol = target * tol_ppm * 1e-6 * 2;
      for (double m : fragment_masses)
      {
        if (m > target / 2 + tol) break;
        const double t = target - m;
        if (t - m < 10.0) continue;
        auto it = std::lower_bound(fragment_masses.begin(), fragment_masses.end(), t - tol);
        if (it != fragment_masses.end() && *it <= t + tol) n++;
      }
    }
    return n;
  };
  static const std::vector<double> shifts {4.3,  -4.3,  6.7,  -6.7,  9.1,  -9.1,  11.3, -11.3, 13.7, -13.7,
                                           17.9, -17.9, 21.3, -21.3, 25.1, -25.1, 29.7, -29.7, 33.3, -33.3};
  std::vector<double> score(cand_masses.size(), 0);
  for (Size i = 0; i < cand_masses.size(); i++)
  {
    const double mass = cand_masses[i];
    double chance = 0;
    for (double s : shifts) chance += pairs_at(mass + s);
    score[i] = pairs_at(mass) - chance / (double)shifts.size();
  }
  // the default is the registered precursor among the target masses (the candidate of the highest charge SNR, chosen without fragments
  // by scan_survey in findPrecursorPeakGroupsForMSnSpectra_; the last such candidate on an exact SNR tie), and the chosen precursor
  // against the decoy masses; another mass replaces it only with at least two pairs beyond chance and two more than the default has
  Size def = cand_masses.size();
  for (Size i = 0; i < cand_masses.size(); i++)
  {
    if (std::abs(cand_masses[i] - default_mass) < 1e-6)
    {
      def = i;
      break;
    }
  }
  if (def == cand_masses.size()) return -1;
  int best = (int)def;
  for (Size i = 0; i < cand_masses.size(); i++)
  {
    if (i == def) continue;
    if (score[i] >= 2.0 && score[i] >= score[def] + 2.0 && score[i] > score[best] + 1e-9) best = (int)i;
  }
  return best == (int)def ? -1 : best;
}

void FLASHDeconvAlgorithm::findPrecursorPeakGroupsForMSnSpectra_(const MSExperiment& map,
                                                                 const std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                                 uint ms_level)
{
  auto cv_of = [](const MSSpectrum& sp) {
    auto fs = sp.getMetaValue("filter string").toString();
    Size p = fs.find("cv=");
    if (p == std::string::npos) return 1e5;
    Size e = fs.find(" ", p);
    if (e == std::string::npos) e = fs.length() - 1;
    return std::stod(fs.substr(p + 3, e - p));
  };
  // the CV of every spectrum of this level and the one below, parsed once (the search below reads them for many neighbours of each scan)
  std::vector<double> spectrum_cvs(map.size(), 0.0);
  std::vector<std::exception_ptr> cv_errors(map.size());
#pragma omp parallel for schedule(dynamic, 256)
  for (SignedSize i = 0; i < (SignedSize)map.size(); i++)
  {
    if (map[i].getMSLevel() != ms_level && map[i].getMSLevel() != ms_level - 1) continue;
    try { spectrum_cvs[i] = cv_of(map[i]); }
    catch (...) { cv_errors[i] = std::current_exception(); }
  }
  for (Size index = 0; index < map.size(); index++)
  {
    const auto& spec = map[index]; // MS2 index
    if (spec.getMSLevel() != ms_level) { continue; }

    int scan_number = getScanNumber(map, index);
    std::string native_id = spec.getNativeID();

    // ADR-0046: a commanded MS2 is located by its command. Null for an uncommanded scan, for MSn > 2, and when
    // the survey the command names is not in this data file (an RT crop) -- those keep the search below unchanged.
    const ScanCommandJoin::Located* located = nullptr;
    if (ms_level == 2)
    {
      auto lit = located_by_scan_.find(scan_number);
      if (lit != located_by_scan_.end() && lit->second.parent_scan_number > 0) { located = &lit->second; }
    }

    // as the serial parse would: the first spectrum (in this order) whose CV does not parse throws
    const auto cv_at = [&](Size i) {
      if (cv_errors[i]) std::rethrow_exception(cv_errors[i]);
      return spectrum_cvs[i];
    };
    const double cv = cv_at(index);

    // The survey (MS n-1) spectra searched for the precursor. A located scan reads the survey its COMMAND names, however many newer
    // surveys were acquired in between (ADR-0046; 88 % of commanded MS2 have at least one), and only that one. Otherwise the preceding
    // survey of the same CV, then the neighbours by distance (one before, one after, two before, ...) up to precursor_rescue_window_
    // away: an LC-MS precursor elutes over several survey scans, so a mass missing from the preceding scan is often in a neighbour
    // (TopFD takes its precursors from LC-MS features).
    std::vector<Size> survey_indices;
    if (located != nullptr)
    {
      for (int i = (int)index - 1; i >= 0; i--)
      {
        if (map[i].getMSLevel() == ms_level - 1 && getScanNumber(map, (Size)i) == located->parent_scan_number)
        {
          survey_indices.push_back((Size)i);
          break;
        }
      }
    }
    else
    {
      std::vector<Size> before, after;
      for (int i = (int)index - 1; i >= 0 && (int)before.size() <= precursor_rescue_window_; i--)
      {
        if (map[i].getMSLevel() == ms_level - 1 && cv_at((Size)i) == cv) before.push_back((Size)i);
      }
      for (Size i = index + 1; i < map.size() && (int)after.size() < precursor_rescue_window_; i++)
      {
        if (map[i].getMSLevel() == ms_level - 1 && cv_at(i) == cv) after.push_back(i);
      }
      if (! before.empty())
      {
        survey_indices.push_back(before[0]);
        for (int k = 1; k <= precursor_rescue_window_; k++)
        {
          if ((Size)k < before.size()) survey_indices.push_back(before[k]);
          if ((Size)(k - 1) < after.size()) survey_indices.push_back(after[k - 1]);
        }
      }
    }
    if (survey_indices.empty()) continue;

    // the deconvolved spectrum of a survey spectrum, by its RT (deconvolved_spectra are sorted by scan number)
    auto deconvolved_of = [&](Size survey_index) -> const DeconvolvedSpectrum* {
      const double rt = map[survey_index].getRT();
      auto it = std::lower_bound(deconvolved_spectra.begin(), deconvolved_spectra.end(), DeconvolvedSpectrum(getScanNumber(map, survey_index) + 1));
      while (it != deconvolved_spectra.begin())
      {
        --it;
        if (std::abs(it->getOriginalSpectrum().getRT() - rt) < 1e-5) return &*it;
        if (it->getOriginalSpectrum().getRT() < rt - 1e-5 && it->getOriginalSpectrum().getMSLevel() == ms_level - 1) break;
      }
      return nullptr;
    };

    // the isolation window
    double start_mz = 0;
    double end_mz = 0;
    double window_target_mz = 0;
    for (auto& precursor : spec.getPrecursors())
    {
      double loffset = precursor.getIsolationWindowLowerOffset();
      double uoffest = precursor.getIsolationWindowUpperOffset();
      loffset = loffset <= 0 ? isolation_window_size_ / 2.0 : loffset;
      uoffest = uoffest <= 0 ? isolation_window_size_ / 2.0 : uoffest;

      // mzML isolation window offsets are relative to the isolation window target m/z. MzMLHandler stores the selected ion m/z as the
      // precursor m/z and, when the two differ, keeps the target as a meta value.
      double target_mz = precursor.metaValueExists("isolation window target m/z") ? (double)precursor.getMetaValue("isolation window target m/z")
                                                                                  : precursor.getMZ();
      window_target_mz = target_mz;
      start_mz = loffset > 100.0 ? loffset : -loffset + target_mz;
      end_mz = uoffest > 100.0 ? uoffest : uoffest + target_mz;
    }

    // The best precursor of one survey spectrum: among the target masses whose peaks lie in the window, the one of the highest charge
    // SNR (at the charge the window implies); for a located scan the COMMANDED mass wins when the survey holds it (ADR-0046: same
    // isotope only, so an isotope-off call competes as an ordinary candidate). With commit, it is registered for this MSn spectrum
    // together with the masses in the window and their intensities, the window bounds, and the target masses of the window as the
    // candidates for the choice by the fragments (with precursor_by_fragments). With decoy reporting, the decoy masses of the window are
    // registered as the challengers of the chosen precursor for the decoy precursor (with precursor_by_fragments, see
    // runSpectralDeconvolution_), or else the one of the highest charge SNR as the decoy precursor. Returns the SNR of the best target
    // (0 or more), or -1 if the window holds none.
    // The raw intensity of a mass's isotope comb at charge c in the isolation window of a survey spectrum (the nearest raw peak within
    // 5 ppm of each isotope position), and whether one of its isotope positions is the isolation target peak (a raw peak within 10 ppm of
    // the target m/z): the window evidence of a candidate, read from the raw spectrum rather than from the deconvolved peak group.
    auto raw_comb = [&](const MSSpectrum& raw, double m, int c) {
      double w = 0;
      if (raw.empty() || c <= 0) return w;
      for (int k = 0; k < 80; k++)
      {
        const double x = (m + k * Constants::ISOTOPE_MASSDIFF_55K_U) / c + Constants::PROTON_MASS_U;
        if (x < start_mz) continue;
        if (x > end_mz) break;
        const Size j = raw.findNearest(x);
        if (std::abs(raw[j].getMZ() - x) <= 5e-6 * x) w += raw[j].getIntensity();
      }
      return w;
    };
    auto owns_target = [&](const MSSpectrum& raw, double m, int c) {
      const double t = window_target_mz;
      if (raw.empty() || c <= 0 || t <= 0) return false;
      const Size j = raw.findNearest(t);
      if (std::abs(raw[j].getMZ() - t) > 10e-6 * t) return false;
      const int k = (int)std::round(((t - Constants::PROTON_MASS_U) * c - m) / Constants::ISOTOPE_MASSDIFF_55K_U);
      if (k < 0) return false;
      const double x = (m + k * Constants::ISOTOPE_MASSDIFF_55K_U) / c + Constants::PROTON_MASS_U;
      return std::abs(x - t) <= 10e-6 * t;
    };

    auto scan_survey = [&](Size survey_index, bool commit) -> double {
      const DeconvolvedSpectrum* ds = deconvolved_of(survey_index);
      if (ds == nullptr) return -1.0;
      const MSSpectrum& precursor_raw_spec = map[survey_index];
      double sum_intensity = .0;
      if (! precursor_raw_spec.empty())
      {
        const Size k = precursor_raw_spec.findNearest(start_mz);
        for (Size l = k + 1; l < precursor_raw_spec.size(); l++)
        {
          if (precursor_raw_spec[l].getMZ() < start_mz) continue;
          if (precursor_raw_spec[l].getMZ() > end_mz) break;
          sum_intensity += precursor_raw_spec[l].getIntensity();
        }
        for (Size l = k;; l--)
        {
          if (precursor_raw_spec[l].getMZ() < start_mz) break;
          if (precursor_raw_spec[l].getMZ() <= end_mz) { sum_intensity += precursor_raw_spec[l].getIntensity(); }
          if (l == 0) break;
        }
      }

      double max_snr = -1.0, max_decoy_snr = -1.0;
      bool found_commanded = false; // ADR-0046: the commanded mass outranks every other species in this window
      double best_residual = 0;
      std::map<double, double> mass_to_intensity;
      mass_to_intensity[.0] = sum_intensity; // total intensity
      const PeakGroup* best_pg = nullptr;
      const FLASHHelperClasses::LogMzPeak* best_peak = nullptr;
      double decoy_mass = 0, decoy_qscore = 0, decoy_qscore1d = 0;
      bool decoy_noise = false;
      // the window bounds (the heaviest target mass and the highest charge at which a target mass has a peak in the window), and, only
      // with commit and precursor_by_fragments, the precursor candidates (the target masses in the window) and their challengers for
      // the decoy precursor (the decoy masses in the window with their noise flags; they do not widen the window bounds)
      double window_max_mass = 0;
      int window_max_charge = 0;
      std::vector<std::pair<PeakGroup, Precursor>> window_cands;
      std::vector<std::pair<double, bool>> window_decoys;
      std::vector<double> window_decoy_snrs; // charge SNR of each window decoy (to order them, the highest first)
      std::vector<double> window_decoy_qscores, window_decoy_qscores1d; // and its precursor score (2D, 1D)
      // every target mass of the window with its most intense peak there and the charge the window implies
      std::vector<std::tuple<const PeakGroup*, const FLASHHelperClasses::LogMzPeak*, int>> window_masses;

      for (auto& pg : *ds)
      {
        if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz) { continue; }
        // With decoy reporting the survey spectrum also holds decoy masses (noise and signal decoys), each in the window if it has a
        // peak there at the charge the window implies. With precursor_by_fragments, all of them challenge the precursor the fragments
        // chose, and one that beats it becomes the decoy precursor of the MSn spectrum (a decoy spectrum with that mass instead of the
        // precursor mass, e.g. for FLASHTnT -decoy_precursor); without, the one of the highest charge SNR becomes it unchecked. Only
        // target masses are precursors.
        if (pg.getTargetDecoyType() != PeakGroup::target)
        {
          int c = int(round(pg.getMonoMass() / start_mz));
          const FLASHHelperClasses::LogMzPeak* top = nullptr;
          for (auto& tmp_peak : pg)
          {
            if (tmp_peak.abs_charge != c || tmp_peak.mz < start_mz || tmp_peak.mz > end_mz) { continue; }
            if (top == nullptr || tmp_peak.intensity > top->intensity) top = &tmp_peak;
          }
          if (top == nullptr) { continue; }
          const bool noise = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy;
          const double snr = pg.getChargeSNR(top->abs_charge);
          if (precursor_by_fragments_)
          {
            if (commit)
            {
              window_decoys.emplace_back(pg.getMonoMass(), noise);
              window_decoy_snrs.push_back(snr);
              window_decoy_qscores.push_back(pg.getQscore2D());
              window_decoy_qscores1d.push_back(pg.getQscore());
            }
            continue;
          }
          if (snr < max_decoy_snr) { continue; }
          max_decoy_snr = snr;
          decoy_mass = pg.getMonoMass();
          decoy_noise = noise;
          decoy_qscore = pg.getQscore2D();
          decoy_qscore1d = pg.getQscore();
          continue;
        }

        double max_intensity = .0;
        const FLASHHelperClasses::LogMzPeak* tmp_precursor = nullptr;

        int c = int(round(pg.getMonoMass() / start_mz));
        double intensity = 0;
        for (auto& tmp_peak : pg)
        {
          if (tmp_peak.abs_charge != c) { continue; }
          if (tmp_peak.mz < start_mz || tmp_peak.mz > end_mz) { continue; }
          intensity += tmp_peak.intensity;

          if (tmp_peak.intensity < max_intensity) { continue; }
          max_intensity = tmp_peak.intensity;
          tmp_precursor = &tmp_peak;
        }

        if (tmp_precursor == nullptr) { continue; }
        mass_to_intensity[pg.getMonoMass()] = intensity;
        window_masses.emplace_back(&pg, tmp_precursor, c);
        window_max_mass = std::max(window_max_mass, pg.getMonoMass());
        window_max_charge = std::max(window_max_charge, c);
        if (commit && precursor_by_fragments_)
        {
          Precursor cand_precursor;
          cand_precursor.setCharge(c);
          cand_precursor.setMZ(tmp_precursor->mz);
          cand_precursor.setIntensity(tmp_precursor->intensity);
          window_cands.emplace_back(pg, cand_precursor);
        }

        // ADR-0046: for a located scan the COMMANDED mass wins when this survey holds it -- SAME isotope only, so an
        // isotope-off call no longer stands in for it and competes as an ordinary candidate. When the window holds
        // no commanded mass the loudest species is the better answer, so this falls through to the charge-SNR rule
        // below and peakGroupFromCommand_ is left for a window that holds nothing at all.
        bool commanded = false;
        if (located != nullptr)
        {
          double residual = 0;
          commanded = ScanCommandJoin::massRank(pg.getMonoMass(), located->row.mono_mass, residual) == 0;
          if (commanded)
          {
            if (found_commanded && residual >= best_residual) { continue; }
            found_commanded = true;
            best_residual = residual;
          }
        }
        if (! commanded)
        {
          if (found_commanded) { continue; } // never displace the commanded mass with a louder neighbour
          auto snr = pg.getChargeSNR(tmp_precursor->abs_charge); // the highest snr one should determine the mass
          if (snr < max_snr) { continue; }
          max_snr = snr;
        }
        best_pg = &pg;
        best_peak = tmp_precursor;
      }

      if (! commit) return best_pg == nullptr ? -1.0 : std::max(max_snr, 0.0);
      // The charge SNR is computed over the whole envelope and cannot tell which ion of the window was isolated. Another mass of the
      // window replaces the choice when its raw comb holds at least 3x the intensity and it is either the same species 14-18 Da away
      // (an oxidation or methylation partner on one comb: the charge SNR took the heavier form in 70 of 71 such TopPIC-only cases) or
      // the owner of the isolation target peak that the choice does not own (Eclip_1476 scan 4033: 4026.21 Da z7 owns the target, share
      // 0.85; the charge SNR took 6311.40 Da z11, share 0.06). On the lysates: 4333 -> 4341 proteoforms at 1%, 3861 -> 3887 with a mass;
      // of 411 changed precursors, raw MS1/MS2 support gained 130 and lost 5. Not for a commanded mass.
      if (best_pg != nullptr && ! found_commanded && window_masses.size() > 1)
      {
        const double m0 = best_pg->getMonoMass();
        const int c0 = int(round(m0 / start_mz));
        const double w0 = raw_comb(precursor_raw_spec, m0, c0);
        const std::tuple<const PeakGroup*, const FLASHHelperClasses::LogMzPeak*, int>* top = nullptr;
        double top_w = 0;
        for (const auto& wm : window_masses)
        {
          const double w = raw_comb(precursor_raw_spec, std::get<0>(wm)->getMonoMass(), std::get<2>(wm));
          if (w > top_w)
          {
            top_w = w;
            top = &wm;
          }
        }
        if (top != nullptr && std::get<0>(*top) != best_pg && top_w >= 3 * w0)
        {
          const double m1 = std::get<0>(*top)->getMonoMass();
          const double d = std::abs(m1 - m0);
          if ((d >= 13.5 && d <= 18.5) || (owns_target(precursor_raw_spec, m1, std::get<2>(*top)) && ! owns_target(precursor_raw_spec, m0, c0)))
          {
            best_pg = std::get<0>(*top);
            best_peak = std::get<1>(*top);
            switched_precursors_++;
          }
        }
      }
      if (precursor_by_fragments_)
      {
        // the decoy precursor is decided by the fragments (runSpectralDeconvolution_). A decoy mass within three isotopes of a target mass
        // of the window is that mass with another monoisotopic peak, not a wrong mass: as a decoy precursor it took the fragments of the
        // right proteoform (FLASHExtender corrects isotope errors back), so most decoys that won were such siblings. It is no decoy here.
        // Nor is a charge harmonic of a target mass (1/2, 1/3, 2 or 3 times it, up to three isotopes of the decoy): FLASHExtender fits the
        // right protein to it as a truncated stretch (Eclip_1480 scan 3012: RS14 31-44 at 1590.69, half the precursor 3181.39).
        std::vector<Size> order;
        for (Size k = 0; k < window_decoys.size(); k++)
        {
          bool sibling = false;
          const double d = window_decoys[k].first;
          for (const auto& [m, j] : mass_to_intensity)
          {
            if (m <= 0) continue;
            const double max_diff = std::max(m, d) * tols_[0] * 2e-6;
            for (int iso = -3; iso <= 3 && ! sibling; iso++)
            {
              const double shift = iso * Constants::ISOTOPE_MASSDIFF_55K_U;
              sibling = std::abs(d - m - shift) <= max_diff;
              for (int h = 2; h <= 3 && ! sibling; h++)
              {
                sibling = std::abs(d - m / h - shift) <= max_diff || std::abs(d - m * h - shift) <= max_diff;
              }
            }
            if (sibling) break;
          }
          if (! sibling) order.push_back(k);
        }
        if (order.empty())
        {
          native_id_window_decoys_.erase(native_id);
          native_id_window_decoy_qscores_.erase(native_id);
          native_id_window_decoy_qscores1d_.erase(native_id);
        }
        else
        {
          std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return window_decoy_snrs[a] > window_decoy_snrs[b]; });
          auto& sorted = native_id_window_decoys_[native_id];
          auto& qscores = native_id_window_decoy_qscores_[native_id];
          auto& qscores1d = native_id_window_decoy_qscores1d_[native_id];
          qscores1d.clear();
          sorted.clear();
          qscores.clear();
          for (Size k : order)
          {
            sorted.push_back(window_decoys[k]);
            qscores[window_decoys[k].first] = window_decoy_qscores[k];
            qscores1d[window_decoys[k].first] = window_decoy_qscores1d[k];
          }
        }
      }
      else if (max_decoy_snr >= 0)
      {
        native_id_decoy_precursor_mass_map_[native_id] = decoy_mass;
        native_id_decoy_precursor_noise_map_[native_id] = decoy_noise;
        native_id_decoy_precursor_qscore_map_[native_id] = decoy_qscore;
        native_id_decoy_precursor_qscore1d_map_[native_id] = decoy_qscore1d;
      }
      if (best_pg != nullptr)
      {
        native_id_precursor_peak_group_map_[native_id] = *best_pg;
        native_id_precursor_peak_group_map_[native_id].setWindowBounds(window_max_mass, window_max_charge);
        // only a window of several candidates leaves the fragments a choice
        if (window_cands.size() > 1) { native_id_window_candidates_[native_id] = std::move(window_cands); }
        else { native_id_window_candidates_.erase(native_id); }
        Precursor precursor;
        precursor.setCharge(best_pg->getRepAbsCharge());
        precursor.setMZ(best_peak->mz);
        precursor.setIntensity(best_peak->intensity);
        native_id_precursor_peak_map_[native_id] = precursor;
      }
      auto& mi = native_id_precursor_mass_intensity_map_[native_id];
      mi.clear();
      for (const auto& [m, j] : mass_to_intensity) mi.emplace_back(m, j);
      return best_pg == nullptr ? -1.0 : std::max(max_snr, 0.0);
    };

    // the first survey spectrum as before; the neighbours only when it holds no precursor (uncommanded scans)
    double snr = scan_survey(survey_indices[0], true);
    if (snr < 0 && survey_indices.size() > 1)
    {
      double best_snr = -1.0;
      Size best_index = 0;
      for (Size k = 1; k < survey_indices.size(); k++)
      {
        const double s = scan_survey(survey_indices[k], false);
        if (s > best_snr) // ties keep the nearer scan
        {
          best_snr = s;
          best_index = survey_indices[k];
        }
      }
      if (best_snr >= 0)
      {
        native_id_decoy_precursor_mass_map_.erase(native_id);
        native_id_decoy_precursor_noise_map_.erase(native_id);
        native_id_decoy_precursor_qscore_map_.erase(native_id);
        native_id_decoy_precursor_qscore1d_map_.erase(native_id);
        scan_survey(best_index, true);
        rescued_precursors_++;
      }
    }
    if (located != nullptr && native_id_precursor_peak_group_map_.find(native_id) == native_id_precursor_peak_group_map_.end())
    {
      native_id_precursor_peak_group_map_[native_id] = peakGroupFromCommand_(*located);
    }
  }
  if (switched_precursors_ > 0)
  {
    OPENMS_LOG_INFO << "Precursor masses switched to a window mass with at least 3x the raw intensity (a 14-18 Da partner or the owner of "
                       "the isolation target peak): "
                    << switched_precursors_ << std::endl;
  }
  if (rescued_precursors_ > 0)
  {
    OPENMS_LOG_INFO << "Precursor masses taken from a neighbouring MS" << (ms_level - 1) << " spectrum (none in the preceding one): " << rescued_precursors_
                    << std::endl;
  }
}

void FLASHDeconvAlgorithm::updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level)
{
  // update precursor feature QScores and qvalues. The survey spectra (of the lower MS level) are not changed below: pointers, not copies
  std::map<int, const DeconvolvedSpectrum*> scan_fullscan;

  for (auto& dspec : deconvolved_spectra)
  {
    if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level - 1) continue;
    int scan = dspec.getScanNumber();
    scan_fullscan[scan] = &dspec;
  }

  for (auto& dspec : deconvolved_spectra)
  {
    if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level) continue;
    if (dspec.getPrecursorPeakGroup().empty()) continue;

    auto precursor_pg = dspec.getPrecursorPeakGroup();

    int pscan = precursor_pg.getScanNumber();
    if (scan_fullscan.find(pscan) == scan_fullscan.end()) continue;

    const auto& fullscan = *scan_fullscan[pscan];

    auto iter = std::lower_bound(fullscan.begin(), fullscan.end(), precursor_pg);

    if (iter == fullscan.end()) continue;

    if (precursor_pg.getMonoMass() == iter->getMonoMass())
    {
      precursor_pg.setFeatureIndex(iter->getFeatureIndex());
      precursor_pg.setQscore(iter->getQscore());
      if (iter->getFeatureIndex() > 0) precursor_pg.setQscore2D(iter->getQscore2D());
    }
    else { precursor_pg.setFeatureIndex(0); }
    dspec.setPrecursorPeakGroup(precursor_pg);
    // the decoy precursor comes from the same survey spectrum: its precursor score is the 2D Qscore of that decoy mass there, the score
    // FLASHDeconv compares its target and decoy masses with (Qvalue), so that the decoy spectrum competes with the spectrum on the precursor
    if (dspec.getDecoyPrecursorMass() > 0)
    {
      // (within 1 ppm; when it is not found, the score it had when the window was scanned stays)
      double traced = -1;
      for (const auto& pg : fullscan)
      {
        if (pg.getTargetDecoyType() == PeakGroup::target || std::abs(pg.getMonoMass() - dspec.getDecoyPrecursorMass()) > dspec.getDecoyPrecursorMass() * 1e-6) continue;
        traced = std::max(traced, pg.getQscore2D());
      }
      if (traced >= 0) dspec.setDecoyPrecursorQscore(traced);
    }
  }
}

Param FLASHDeconvAlgorithm::featureTraceParam_() const
{
  Param mf_param = param_.copy("ft:", true);

  if (((double)mf_param.getValue("mass_error_ppm")) < 0) { mf_param.setValue("mass_error_ppm", tols_[0]); }
  if (((double)mf_param.getValue("min_cos")) < 0) { mf_param.setValue("min_cos", min_cos_[0]); }

  mf_param.setValue("noise_threshold_int", .0);
  mf_param.setValue("reestimate_mt_sd", "false");
  mf_param.setValue("trace_termination_criterion", "outlier");
  mf_param.setValue("trace_termination_outliers", 20);
  mf_param.setValue("chrom_peak_snr", .0);
  return mf_param;
}

std::vector<FLASHHelperClasses::MassFeature> FLASHDeconvAlgorithm::traceMinLevelFeatures_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
{
  auto mass_tracer = MassFeatureTrace();
  mass_tracer.setParameters(featureTraceParam_()); // maybe go to set param
  // Find features for MS1 or the minimum MS level in the dataset.
  auto deconvolved_features = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), deconvolved_spectra, (int)current_min_ms_level_, false);

  if (report_decoy_)
  {
    const auto& decoy_deconvolved_features = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), deconvolved_spectra, (int)current_min_ms_level_, true);
    deconvolved_features.insert(deconvolved_features.end(), decoy_deconvolved_features.begin(), decoy_deconvolved_features.end());
  }
  return deconvolved_features;
}

void FLASHDeconvAlgorithm::runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                              std::vector<FLASHHelperClasses::MassFeature>& deconvolved_features)
{
  // the features of the minimum MS level: traced during the wave of the last MS level if possible (runSpectralDeconvolution_)
  if (min_level_trace_error_)
  {
    auto error = min_level_trace_error_;
    min_level_trace_error_ = nullptr;
    std::rethrow_exception(error);
  }
  if (min_level_features_)
  {
    deconvolved_features = std::move(*min_level_features_);
    min_level_features_.reset();
  }
  else { deconvolved_features = traceMinLevelFeatures_(deconvolved_spectra); }

  auto mass_tracer = MassFeatureTrace();
  Param mf_param = featureTraceParam_();
  mf_param.setValue("min_trace_length", 1e-5); // allow all traces for MSn
  mass_tracer.setParameters(mf_param);
  // Find features for MSn
  for (int ms_level = (int)current_min_ms_level_ + 1; ms_level <= (int)current_max_ms_level_; ms_level++)
  {
    updatePrecursorQScores_(deconvolved_spectra, ms_level);
    // An MSn precursor that belongs to a feature of the lowest MS level gets the feature's mass, at the isotope label the precursor has
    // (chosen with the fragments: precursor_by_fragments). The feature mass is the centroid of the masses of all its survey spectra, so it
    // scatters less than the mass of the one survey spectrum the precursor was found in; TopFD likewise gives an MS2 spectrum its
    // feature's mass. Three E. coli lysates, 1% PrSMs whose precursor mass equals the sequence mass: 2.5-97.5 % of the error -0.96..+2.23
    // ppm with the survey mass, -0.48..+1.71 ppm with the feature mass (TopFD -0.41..+1.45); identifications unchanged (5106 -> 5108
    // proteoforms, 1008 -> 1009 proteins). Taking the feature's isotope label as well put one third fewer precursors one isotope low, but
    // left the reported proteoform masses as they were and lost 11 proteoforms.
    if (ms_level == (int)current_min_ms_level_ + 1)
    {
      std::map<uint, double> feature_mass;
      for (const auto& mf : deconvolved_features)
      {
        if (! mf.is_decoy) feature_mass[mf.index] = mf.mt.getCentroidMZ() + mf.iso_offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      }
      for (auto& dspec : deconvolved_spectra)
      {
        if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level || dspec.getPrecursorPeakGroup().empty()) continue;
        auto pg = dspec.getPrecursorPeakGroup();
        if (pg.getFeatureIndex() == 0) continue;
        const auto it = feature_mass.find(pg.getFeatureIndex());
        if (it == feature_mass.end()) continue;
        const double d = pg.getIsotopeDaDistance();
        pg.setMonoisotopicMass(it->second + std::round((pg.getMonoMass() - it->second) / d) * d);
        dspec.setPrecursorPeakGroup(pg);
      }
    }
    std::map<uint, std::vector<Size>> feature_index_set;
    for (Size i = 0; i < deconvolved_spectra.size(); i++)
    {
      const auto& dspec = deconvolved_spectra[i];
      if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level) continue;
      if (dspec.getPrecursorPeakGroup().empty()) continue;
      uint findex = dspec.getPrecursorPeakGroup().getFeatureIndex();
      if (findex == 0) continue;
      feature_index_set[findex].push_back(i);
    }

    // The groups are traced in parallel, each from its own spectra (every spectrum is in one group; the tracing reads no feature index or
    // 2D Qscore), and numbered and labelled afterwards in the group order (the target, then the decoy features of a group), the order in
    // which they were traced and numbered one after the other before.
    std::vector<std::pair<uint, std::vector<Size>>> groups(feature_index_set.begin(), feature_index_set.end());
    std::vector<MassFeatureTrace::TracedFeatures> traced(groups.size()), traced_decoy(report_decoy_ ? groups.size() : 0);
    std::vector<std::exception_ptr> errors(groups.size());
    bool parallel = groups.size() > 1;
#ifdef _OPENMP
    parallel = parallel && ! omp_in_parallel();
#endif
#pragma omp parallel if (parallel) default(none) shared(groups, traced, traced_decoy, errors, deconvolved_spectra, mass_tracer, ms_level)
    {
      MassFeatureTrace tracer(mass_tracer);
#pragma omp for schedule(dynamic, 1)
      for (SignedSize g = 0; g < (SignedSize)groups.size(); g++)
      {
        const auto& indices = groups[g].second;
        std::vector<DeconvolvedSpectrum> tmp_dspec;
        // moved out and back: each spectrum is in one feature group, and nothing else reads it in between
        try
        {
          tmp_dspec.reserve(indices.size());
          for (Size i : indices)
            tmp_dspec.push_back(std::move(deconvolved_spectra[i]));
          traced[g] = tracer.traceFeatures(sd_.getAveragine(), tmp_dspec, ms_level, false);
          if (report_decoy_) traced_decoy[g] = tracer.traceFeatures(sd_.getAveragine(), tmp_dspec, ms_level, true);
        }
        catch (...)
        {
          errors[g] = std::current_exception();
        }
        for (Size j = 0; j < tmp_dspec.size(); j++)
          deconvolved_spectra[indices[j]] = std::move(tmp_dspec[j]);
      }
    }
    for (Size g = 0; g < groups.size(); g++)
    {
      if (errors[g]) std::rethrow_exception(errors[g]);
      auto df = MassFeatureTrace::numberFeatures(traced[g], deconvolved_spectra, groups[g].second);
      deconvolved_features.insert(deconvolved_features.end(), std::make_move_iterator(df.begin()), std::make_move_iterator(df.end()));
      if (report_decoy_)
      {
        auto df_decoy = MassFeatureTrace::numberFeatures(traced_decoy[g], deconvolved_spectra, groups[g].second);
        deconvolved_features.insert(deconvolved_features.end(), std::make_move_iterator(df_decoy.begin()), std::make_move_iterator(df_decoy.end()));
      }
    }
  }
}
} // namespace OpenMS