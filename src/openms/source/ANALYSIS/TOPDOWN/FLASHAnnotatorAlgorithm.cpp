// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <unordered_map>

namespace OpenMS
{
namespace
{
  /// internal (residue) monoisotopic masses by one-letter code; 0 for codes without a residue of positive mass
  const std::vector<double>& residueMasses()
  {
    static const std::vector<double> masses = [] {
      std::vector<double> m(256, 0.0);
      const ResidueDB* db = ResidueDB::getInstance();
      for (int c = 'A'; c <= 'Z'; c++)
      {
        const Residue* r = db->getResidue((unsigned char)c);
        if (r == nullptr) continue;
        double w = r->getMonoWeight(Residue::Internal);
        if (w <= 0) continue;
        m[c] = w;
        m[c - 'A' + 'a'] = w;
      }
      return m;
    }();
    return masses;
  }

  // the index of the first peak at or above low in the whole spectrum (the size if none)
  inline int firstPeakAtOrAbove(const FLASHAnnotatorAlgorithm::PeakArrays& spec, double low)
  {
    const double* base = spec.mzs.data();
    const int size = (int)spec.mzs.size();
    const auto& buckets = spec.bucket_start;
    if (buckets.empty()) return (int)(std::lower_bound(base, base + size, low) - base);
    const double b = std::floor(low / FLASHAnnotatorAlgorithm::PeakArrays::bucket_width);
    if (b < 0) return 0;
    if (b + 1 >= (double)buckets.size()) return size;
    return (int)(std::lower_bound(base + buckets[(Size)b], base + buckets[(Size)b + 1], low) - base);
  }

  /**
    Index of the observed peak nearest to @p target within [target - w, target + w], or -1 (the first of equally near
    peaks wins). The first index with m/z >= target - w is searched in [lo, hi), which must contain it (or hi if it is
    hi), and returned in @p lb so that the next isotope peak can be searched in a narrower range.
  */
  // the index of the peak nearest to target within w among peaks lo .. hi - 1, -1 if none; lb is the first index in lo .. hi at or above
  // target - w. With the bucket index of the spectrum, that index is searched in one bucket only: the first peak at or above target - w
  // of the whole spectrum lies between the starts of its bucket and the next, and clamped to lo .. hi it is the one of the range.
  inline int nearestPeak(const FLASHAnnotatorAlgorithm::PeakArrays& spec, int lo, int hi, double target, double w, int& lb)
  {
    const double* base = spec.mzs.data();
    const int size = (int)spec.mzs.size();
    lb = std::max(lo, std::min(hi, firstPeakAtOrAbove(spec, target - w)));
    int best = -1;
    double best_d = 0;
    for (int i = lb; i < size && base[i] <= target + w; i++)
    {
      double d = std::abs(base[i] - target);
      if (best < 0 || d < best_d)
      {
        best = i;
        best_d = d;
      }
    }
    return best;
  }
} // namespace

FLASHAnnotatorAlgorithm::FLASHAnnotatorAlgorithm(): DefaultParamHandler("FLASHAnnotatorAlgorithm")
{
  setDefaultParams_();
}

void FLASHAnnotatorAlgorithm::setDefaultParams_()
{
  defaults_.setValue("tol", 10.0, "Tolerance (ppm) for matching isotope peaks to observed peaks, after calibration.");
  defaults_.setMinFloat("tol", 0.0);

  defaults_.setValue("calibrate", "true",
                     "Correct the observed m/z by the median ppm error of the fragment envelopes matched at calibration_tol "
                     "(computed per PrSM, from its b/y fragments).");
  defaults_.setValidStrings("calibrate", {"true", "false"});

  defaults_.setValue("calibration_tol", 20.0, "Tolerance (ppm) of the calibration pass.");
  defaults_.setMinFloat("calibration_tol", 0.0);

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Ion types to annotate.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaults_.setValue("min_peaks", 3, "Minimum number of matched isotope peaks of an envelope, the most abundant one included.");
  defaults_.setMinInt("min_peaks", 1);

  defaults_.setValue("min_rel_isotope_intensity", 0.1,
                     "Isotope peaks of at least this intensity relative to the most abundant one are predicted.");
  defaults_.setMinFloat("min_rel_isotope_intensity", 0.0);
  defaults_.setMaxFloat("min_rel_isotope_intensity", 1.0);

  defaults_.setValue("min_cos", 0.8, "Minimum cosine between predicted and observed isotope intensities (unmatched peaks count as 0).");
  defaults_.setMinFloat("min_cos", -1.0);
  defaults_.setMaxFloat("min_cos", 1.0);

  defaults_.setValue("min_pearson", 0.7,
                     "Minimum Pearson correlation between predicted and observed isotope intensities (unmatched peaks count as "
                     "0). It rejects flat runs of peaks that pass the cosine.");
  defaults_.setMinFloat("min_pearson", -1.0);
  defaults_.setMaxFloat("min_pearson", 1.0);

  defaults_.setValue("min_mass_per_charge", 300.0, "Fragment charges are capped at floor(fragment mass / this value).");
  defaults_.setMinFloat("min_mass_per_charge", 1.0);

  defaults_.setValue("max_mass", 200000.0, "Fragments heavier than this are not annotated.");
  defaults_.setMinFloat("max_mass", 100.0);

  // whole multiples (9, 13, 21, 33, 39, 53) of 1.000495 Da, the spacing of the peptide mass bands
  defaults_.setValue("decoy_shifts",
                     std::vector<double> {-53.026235, -39.019305, -33.016335, -21.010395, -13.006435, -9.004455, 9.004455, 13.006435,
                                          21.010395, 33.016335, 39.019305, 53.026235},
                     "Mass shifts (Da) of the decoy fragment sets that measure the chance match rate of each PrSM in its own "
                     "spectrum. Whole multiples of 1.000495 Da keep the decoy masses in the peptide mass bands, where random "
                     "fragments lie; shifts between the bands match less often and underestimate the chance rate. Empty: no "
                     "chance rate, no RawScore and no RawLLR.");

  defaults_.setValue("p1", 0.3, "Probability that a fragment of a correct proteoform is observed, for the log-likelihood ratio.");
  defaults_.setMinFloat("p1", 1e-6);
  defaults_.setMaxFloat("p1", 1.0 - 1e-6);

  defaults_.setValue("report_decoy_envelopes", "false",
                     "Also report the accepted envelopes of the decoy fragment sets (every charge), e.g. to study how chance envelopes "
                     "differ from real ones.");
  defaults_.setValidStrings("report_decoy_envelopes", {"true", "false"});
  defaults_.addTag("report_decoy_envelopes", "advanced");

  defaults_.setValue("min_peaks_small", 2,
                     "Fragments with fewer than min_peaks predicted isotope peaks (small fragments) are annotated if all their predicted "
                     "peaks, at least this many, are matched (cosine only: two points always correlate). 0: not annotated.");
  defaults_.setMinInt("min_peaks_small", 0);

  defaults_.setValue("terminal_calibration", std::vector<double> {0.479, 0.697},
                     "c0, c1: -log10 P(RawScore >= s) = c0 + c1 s (s >= 1; linear from 0 below), the tail of RawScore for shuffled "
                     "proteoforms (composition-preserving shuffles of E. coli lysate PrSMs). Used for FAScore.");
  defaults_.setValue("internal_calibration", std::vector<double> {-11.531, 2.062, -0.471},
                     "a, b, c: logit P(K >= k | n) = a + b log10(n) + c k, the chance of k counted internal pieces among n tested ones, "
                     "fitted on shuffled proteoforms. Used for FAScore.");

  defaults_.setValue("internal_ions", "true",
                     "Annotate internal fragments (both ends inside the proteoform). The cleavage after a residue is authentic if pieces of "
                     "consecutive lengths that all end with that residue, or all start right after it, are accepted (a ladder giving a tag): "
                     "a tag of at least internal_min_tag_length_asp residues after Asp, internal_min_tag_length_pro before Pro, "
                     "internal_min_tag_length elsewhere. Ladder pieces "
                     "count; internal_extension decides which other pieces at authentic cleavages count. The same decoy shifts measure "
                     "their chance rate.");
  defaults_.setValidStrings("internal_ions", {"true", "false"});

  defaults_.setValue("internal_max_length", 150, "Longest internal piece (residues) that is tested.");
  defaults_.setMinInt("internal_max_length", 2);

  defaults_.setValue("internal_min_tag_length_asp", 2,
                     "Tag length (consecutive pieces - 1) that makes a cleavage after Asp authentic (aspartic acid effect). On raw peaks "
                     "of the 20260715 targeted set, 1 adds chance coverage of ~2.4 points, 2 ~0.7.");
  defaults_.setMinInt("internal_min_tag_length_asp", 1);

  defaults_.setValue("internal_min_tag_length_pro", 4,
                     "Tag length (consecutive pieces - 1) that makes a cleavage before Pro authentic (proline effect). A cleavage after "
                     "Asp and before Pro takes the smaller of the two thresholds. The value of internal_min_tag_length treats Pro like "
                     "any other residue. On native spectra of 34 E. coli proteins, tags of length 4 before Pro are 4% chance, length 3 9%.");
  defaults_.setMinInt("internal_min_tag_length_pro", 1);

  defaults_.setValue("internal_min_tag_length", 5,
                     "Tag length that makes any other cleavage authentic. On raw peaks of the 20260715 targeted set, 3 adds chance "
                     "coverage of ~1.1 points, 5 ~0.02.");
  defaults_.setMinInt("internal_min_tag_length", 1);

  defaults_.setValue("internal_extension", "tag1",
                     "With internal_asp, which accepted internal pieces with an end at an authentic Asp cleavage count beyond its ladders: "
                     "'none', 'tag1' (pieces with a neighbor one residue longer or shorter sharing the Asp end, i.e. a tag of length "
                     "one) or 'all'. For a PrSM whose start or end position is undetermined (-1), 'none' is used.");
  defaults_.setValidStrings("internal_extension", {"none", "tag1", "all"});

  defaultsToParam_();
}

void FLASHAnnotatorAlgorithm::updateMembers_()
{
  tol_ = param_.getValue("tol");
  calibrate_on_ = param_.getValue("calibrate").toString() == "true";
  cal_tol_ = param_.getValue("calibration_tol");
  min_peaks_ = param_.getValue("min_peaks");
  min_rel_iso_ = param_.getValue("min_rel_isotope_intensity");
  min_cos_ = param_.getValue("min_cos");
  min_pearson_ = param_.getValue("min_pearson");
  min_mass_per_charge_ = param_.getValue("min_mass_per_charge");
  max_mass_ = param_.getValue("max_mass");
  decoy_shifts_ = param_.getValue("decoy_shifts").toDoubleVector();
  p1_ = param_.getValue("p1");
  report_decoys_ = param_.getValue("report_decoy_envelopes").toString() == "true";
  internal_on_ = param_.getValue("internal_ions").toString() == "true";
  min_peaks_small_ = param_.getValue("min_peaks_small");
  terminal_calibration_ = param_.getValue("terminal_calibration").toDoubleVector();
  internal_calibration_ = param_.getValue("internal_calibration").toDoubleVector();
  if (terminal_calibration_.size() != 2) terminal_calibration_ = {0.479, 0.697};
  if (internal_calibration_.size() != 3) internal_calibration_ = {-11.531, 2.062, -0.471};
  internal_max_length_ = param_.getValue("internal_max_length");
  min_run_asp_ = (int)param_.getValue("internal_min_tag_length_asp") + 1;
  min_run_other_ = (int)param_.getValue("internal_min_tag_length") + 1;
  min_run_pro_ = (int)param_.getValue("internal_min_tag_length_pro") + 1;
  const std::string extension = param_.getValue("internal_extension").toString();
  extension_min_run_ = extension == "all" ? 1 : (extension == "tag1" ? 2 : 0);

  // same fragment masses as FLASHTaggerAlgorithm and FLASHExtenderAlgorithm: internal residue sum + ion type shift
  ion_types_.clear();
  ion_shifts_.clear();
  ion_is_prefix_.clear();
  for (const auto& ion_str : param_.getValue("ion_type").toStringVector())
  {
    if (std::find(ion_types_.begin(), ion_types_.end(), ion_str) != ion_types_.end()) continue;
    double shift;
    bool prefix;
    if (ion_str == "a") { shift = Residue::getInternalToAIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "b") { shift = Residue::getInternalToBIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "c") { shift = Residue::getInternalToCIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "x") { shift = Residue::getInternalToXIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "y") { shift = Residue::getInternalToYIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "z") { shift = Residue::getInternalToZIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "zp1") { shift = Residue::getInternalToZp1Ion().getMonoWeight(); prefix = false; }
    else if (ion_str == "zp2") { shift = Residue::getInternalToZp2Ion().getMonoWeight(); prefix = false; }
    else { continue; }
    ion_types_.push_back(ion_str);
    ion_shifts_.push_back(shift);
    ion_is_prefix_.push_back(prefix);
  }

  // calibration uses the b and y fragments (as the prototype did); all ion types if neither is given
  bool has_by = std::find(ion_types_.begin(), ion_types_.end(), "b") != ion_types_.end()
                || std::find(ion_types_.begin(), ion_types_.end(), "y") != ion_types_.end();
  ion_for_calibration_.clear();
  for (const auto& ion : ion_types_)
  {
    ion_for_calibration_.push_back(! has_by || ion == "b" || ion == "y");
  }

  // ion types of the same terminus whose shifts differ by a few isotope spacings (z, z+1, z+2) match the same envelope one
  // or two isotope peaks apart: they form one fragment hypothesis for the score
  const int ion_count = (int)ion_types_.size();
  ion_group_.assign(ion_count, -1);
  ion_group_count_ = 0;
  for (int t = 0; t < ion_count; t++)
  {
    if (ion_group_[t] >= 0) continue;
    ion_group_[t] = ion_group_count_;
    // grow the group of t: an ion type joins it if it is equivalent to any member
    for (bool grown = true; grown;)
    {
      grown = false;
      for (int u = t + 1; u < ion_count; u++)
      {
        if (ion_group_[u] >= 0 || ion_is_prefix_[u] != ion_is_prefix_[t]) continue;
        for (int v = t; v < ion_count; v++)
        {
          if (ion_group_[v] != ion_group_count_) continue;
          const double diff = ion_shifts_[u] - ion_shifts_[v];
          const int offset = (int)std::round(diff / Constants::C13C12_MASSDIFF_U);
          if (offset == 0 || std::abs(offset) > max_isotope_offset_) continue;
          if (std::abs(diff - offset * Constants::C13C12_MASSDIFF_U) > isotope_equivalence_tol_) continue;
          ion_group_[u] = ion_group_count_;
          grown = true;
          break;
        }
      }
    }
    ion_group_count_++;
  }

  // internal ion types: an N-end from a suffix-type cleavage and a C-end from a prefix-type cleavage (b/y: residue sum)
  internal_offsets_.clear();
  internal_ion_types_.clear();
  const double water = Residue::getInternalToFull().getMonoWeight();
  for (int c = 0; c < ion_count; c++)
  {
    if (! ion_is_prefix_[c]) continue;
    for (int n = 0; n < ion_count; n++)
    {
      if (ion_is_prefix_[n]) continue;
      const double off = ion_shifts_[c] + ion_shifts_[n] - water;
      if (std::any_of(internal_offsets_.begin(), internal_offsets_.end(), [off](double o) { return std::abs(o - off) < 1e-4; })) continue;
      internal_offsets_.push_back(off);
      internal_ion_types_.push_back(ion_types_[c] + "/" + ion_types_[n]);
    }
  }
}

const std::vector<std::string>& FLASHAnnotatorAlgorithm::getIonTypes() const
{
  return ion_types_;
}

const std::vector<std::string>& FLASHAnnotatorAlgorithm::getInternalIonTypes() const
{
  return internal_ion_types_;
}

int FLASHAnnotatorAlgorithm::requiredPeaks_(int predicted) const
{
  if (predicted >= min_peaks_) return min_peaks_;
  // small fragments, with fewer predicted isotope peaks: all of them
  if (min_peaks_small_ > 0 && predicted >= min_peaks_small_) return predicted;
  return -1;
}

double FLASHAnnotatorAlgorithm::calibratedTerminalLog10P_(double score) const
{
  // the shuffle-null tail of RawScore: -log10 P(RawScore >= s) = c0 + c1 s above s = 1, linear from 0 below
  if (score <= 0) return 0;
  const double at_one = terminal_calibration_[0] + terminal_calibration_[1];
  return score < 1 ? score * at_one : terminal_calibration_[0] + terminal_calibration_[1] * score;
}

double FLASHAnnotatorAlgorithm::internalLog10P_(int k, int n) const
{
  if (k <= 0 || n <= 0) return 0;
  // logit P(K >= k | n) = a + b log10(n) + c k, from internal pieces counted for shuffled proteoforms
  const double z = internal_calibration_[0] + internal_calibration_[1] * std::log10((double)n) + internal_calibration_[2] * k;
  // -log10 of 1 / (1 + exp(-z)), computed without overflow
  const double log_p = z >= 0 ? -std::log1p(std::exp(-z)) : z - std::log1p(std::exp(z));
  return std::max(0.0, -log_p / std::log(10.0));
}

double FLASHAnnotatorAlgorithm::densityAt_(const PeakArrays& spec, double mz)
{
  auto lo = std::lower_bound(spec.mzs.begin(), spec.mzs.end(), mz - 1.0);
  auto hi = std::upper_bound(spec.mzs.begin(), spec.mzs.end(), mz + 1.0);
  return (double)(hi - lo) / 2.0;
}

double FLASHAnnotatorAlgorithm::rankOf_(const std::vector<double>& sorted_intensities, double intensity)
{
  if (sorted_intensities.empty()) return 0.0;
  return (double)(std::upper_bound(sorted_intensities.begin(), sorted_intensities.end(), intensity) - sorted_intensities.begin())
         / (double)sorted_intensities.size();
}

void FLASHAnnotatorAlgorithm::PeakArrays::indexBuckets()
{
  bucket_start.clear();
  if (mzs.empty()) return;
  const Size count = (Size)std::floor(mzs.back() / bucket_width) + 2;
  bucket_start.resize(count);
  Size i = 0;
  for (Size b = 0; b < count; b++)
  {
    while (i < mzs.size() && mzs[i] < (double)b * bucket_width) i++;
    bucket_start[b] = (int)i;
  }
}

void FLASHAnnotatorAlgorithm::PeakOccupancy::build(const PeakArrays& spec)
{
  bits.clear();
  if (spec.mzs.empty()) return;
  first_mz = spec.mzs.front();
  last_mz = spec.mzs.back();
  // bins_per_th is a power of two: mz * bins_per_th is exact, and the bin of a peak in [low, high] lies between the bins of low and high
  first_bin = std::floor(first_mz * bins_per_th);
  const double bins = std::floor(last_mz * bins_per_th) - first_bin + 1;
  // not built for an unusual m/z range; nor with negative m/z, where the bucket search of nearestPeak starts at the first peak
  if (! (first_mz >= 0 && bins >= 1 && bins <= (double)(1 << 24))) return;
  bits.assign(((Size)bins + 63) / 64, 0);
  double previous = first_mz;
  for (const double mz : spec.mzs)
  {
    if (! (mz >= previous && mz <= last_mz)) // not sorted: not built
    {
      bits.clear();
      return;
    }
    previous = mz;
    const Size b = (Size)(std::floor(mz * bins_per_th) - first_bin);
    bits[b >> 6] |= UInt64(1) << (b & 63);
  }
}

bool FLASHAnnotatorAlgorithm::PeakOccupancy::mayHavePeak(double low, double high) const
{
  if (bits.empty()) return true;
  // a peak in [low, high] lies in [l, h], the window within the peaks' m/z range
  const double l = std::max(low, first_mz);
  const double h = std::min(high, last_mz);
  if (l > h) return false;
  if (! (l <= h)) return true; // not a number: searched
  const Size lo = (Size)(std::floor(l * bins_per_th) - first_bin);
  const Size hi = (Size)(std::floor(h * bins_per_th) - first_bin);
  for (Size word = lo >> 6; word <= hi >> 6; word++)
  {
    UInt64 m = bits[word];
    if (word == lo >> 6) m &= ~UInt64(0) << (lo & 63);
    if (word == hi >> 6) m &= ~UInt64(0) >> (63 - (hi & 63));
    if (m != 0) return true;
  }
  return false;
}

void FLASHAnnotatorAlgorithm::toPeakArrays(const MSSpectrum& spec, PeakArrays& peaks)
{
  peaks.mzs.clear();
  peaks.intensities.clear();
  peaks.mzs.reserve(spec.size());
  peaks.intensities.reserve(spec.size());

  if (spec.isSorted())
  {
    for (const auto& p : spec)
    {
      if (p.getIntensity() <= 0) continue;
      peaks.mzs.push_back(p.getMZ());
      peaks.intensities.push_back(p.getIntensity());
    }
    peaks.indexBuckets();
    return;
  }

  std::vector<std::pair<double, double>> sorted;
  sorted.reserve(spec.size());
  for (const auto& p : spec)
  {
    if (p.getIntensity() <= 0) continue;
    sorted.emplace_back(p.getMZ(), p.getIntensity());
  }
  std::sort(sorted.begin(), sorted.end());
  for (const auto& [mz, intensity] : sorted)
  {
    peaks.mzs.push_back(mz);
    peaks.intensities.push_back(intensity);
  }
  peaks.indexBuckets();
}

double FLASHAnnotatorAlgorithm::getProteoformMass(const std::string& sequence, const std::vector<Modification>& modifications)
{
  static const double water = Residue::getInternalToFull().getMonoWeight();
  const auto& aa = residueMasses();
  double mass = water;
  for (const char c : sequence)
  {
    double m = aa[(unsigned char)c];
    if (m <= 0) return -1;
    mass += m;
  }
  for (const auto& mod : modifications)
  {
    mass += mod.mass;
  }
  return mass;
}

void FLASHAnnotatorAlgorithm::ensureIsotopeTable_(double mass)
{
  mass = std::min(mass, max_mass_);
  if (table_min_rel_iso_ == min_rel_iso_ && table_.max_mass >= mass && ! table_.apex_pos.empty()) return;

  // bins at min_mass_ + b * bin_width_ (50, 75, ... Da) as PrecalculatedAveragine(50, max mass, 25) in
  // SpectralDeconvolution::calculateAveragine, up to the bin nearest to mass
  const int bin_count = std::max(2, (int)std::floor((mass - min_mass_) / bin_width_) + 2);
  const double top_mass = min_mass_ + (bin_count - 1) * bin_width_;

  CoarseIsotopePatternGenerator generator(300);
  auto max_iso = generator.estimateFromPeptideWeight(top_mass);
  max_iso.trimRight(0.0001 * max_iso.getMostAbundant().getIntensity());
  const int max_isotope = std::max(200, (int)max_iso.size());

  std::vector<std::vector<int>> iso_index(bin_count);
  std::vector<std::vector<double>> iso_rel(bin_count);
  std::vector<int> apex_pos(bin_count, 0), mono_pos(bin_count, -1);

  // the distribution of a bin does not depend on the other bins, so chunks of bins are computed in parallel, each by its
  // own PrecalculatedAveragine. Each bin is then flattened into its kept peaks.
  const int chunk_size = 100;
  const int chunk_count = (bin_count + chunk_size - 1) / chunk_size;
#pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < chunk_count; c++)
  {
    const int first = c * chunk_size;
    const int last = std::min(bin_count, first + chunk_size) - 1;
    CoarseIsotopePatternGenerator chunk_generator(max_isotope);
    FLASHHelperClasses::PrecalculatedAveragine avg(min_mass_ + first * bin_width_, min_mass_ + last * bin_width_, bin_width_, chunk_generator, false);
    for (int b = first; b <= last; b++)
    {
      const IsotopeDistribution iso = avg.get(min_mass_ + b * bin_width_);
      Size apex = 0;
      for (Size k = 1; k < iso.size(); k++)
      {
        if (iso[k].getIntensity() > iso[apex].getIntensity()) apex = k;
      }
      const double apex_intensity = iso.size() == 0 ? 0.0 : (double)iso[apex].getIntensity();
      if (apex_intensity <= 0) continue;
      for (Size k = 0; k < iso.size(); k++)
      {
        const double rel = iso[k].getIntensity() / apex_intensity;
        if (rel <= 0 || rel < min_rel_iso_) continue;
        if (k == apex) apex_pos[b] = (int)iso_rel[b].size();
        if (k == 0) mono_pos[b] = 0;
        iso_index[b].push_back((int)k);
        iso_rel[b].push_back(rel);
      }
    }
  }

  table_ = IsotopeTable();
  table_.max_mass = top_mass;
  table_.offset.resize(bin_count + 1, 0);
  table_.apex_pos = apex_pos;
  table_.mono_pos = mono_pos;
  table_.rel_sum.resize(bin_count, 0);
  table_.rel_sum_sq.resize(bin_count, 0);
  for (int b = 0; b < bin_count; b++)
  {
    table_.offset[b] = (int)table_.iso_index.size();
    for (Size j = 0; j < iso_rel[b].size(); j++)
    {
      table_.iso_index.push_back(iso_index[b][j]);
      table_.iso_rel.push_back(iso_rel[b][j]);
      table_.rel_sum[b] += iso_rel[b][j];
      table_.rel_sum_sq[b] += iso_rel[b][j] * iso_rel[b][j];
    }
  }
  table_.offset[bin_count] = (int)table_.iso_index.size();
  table_min_rel_iso_ = min_rel_iso_;
}

int FLASHAnnotatorAlgorithm::massToBin_(double mass) const
{
  // the same rounding as PrecalculatedAveragine::massToIndex_
  int bin = (int)std::round(std::max(.0, mass - min_mass_) / bin_width_);
  return std::min(bin, (int)table_.apex_pos.size() - 1);
}

void FLASHAnnotatorAlgorithm::getIsotopeEnvelope(double mono_mass, int charge, std::vector<double>& mzs, std::vector<double>& intensities)
{
  mzs.clear();
  intensities.clear();
  if (charge <= 0 || mono_mass <= 0) return;
  ensureIsotopeTable_(mono_mass + bin_width_);
  const int bin = massToBin_(mono_mass);
  for (int j = table_.offset[bin]; j < table_.offset[bin + 1]; j++)
  {
    mzs.push_back((mono_mass + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U + charge * Constants::PROTON_MASS_U) / charge);
    intensities.push_back(table_.iso_rel[j]);
  }
}

bool FLASHAnnotatorAlgorithm::getFragments(const PrSM& prsm, std::vector<Fragment>& fragments, std::vector<Fragment>* span_fragments) const
{
  fragments.clear();
  if (span_fragments != nullptr) span_fragments->clear();
  const auto& aa = residueMasses();
  const int length = (int)prsm.sequence.size();
  std::vector<double> prefix_mass(length + 1, .0);
  for (int i = 0; i < length; i++)
  {
    double m = aa[(unsigned char)prsm.sequence[i]];
    if (m <= 0) return false;
    prefix_mass[i + 1] = prefix_mass[i] + m;
  }

  // modification ranges, clamped into the proteoform
  std::vector<std::pair<int, int>> ranges;
  ranges.reserve(prsm.modifications.size());
  for (const auto& mod : prsm.modifications)
  {
    int s = std::max(0, std::min(length - 1, std::min(mod.start, mod.end)));
    int e = std::max(0, std::min(length - 1, std::max(mod.start, mod.end)));
    ranges.emplace_back(s, e);
  }

  fragments.reserve(2 * ion_types_.size() * (Size)std::max(0, length));
  for (int t = 0; t < (int)ion_types_.size(); t++)
  {
    const bool prefix = ion_is_prefix_[t];
    if ((prefix && ! prsm.known_n_term) || (! prefix && ! prsm.known_c_term)) continue;
    for (int len = 1; len < length; len++)
    {
      // residues [first, last] of the fragment
      const int first = prefix ? 0 : length - len;
      const int last = prefix ? len - 1 : length - 1;
      double mod_mass = 0;
      bool split = false;
      for (Size m = 0; m < ranges.size(); m++)
      {
        if (ranges[m].first >= first && ranges[m].second <= last) { mod_mass += prsm.modifications[m].mass; }
        else if (ranges[m].second >= first && ranges[m].first <= last)
        {
          split = true; // the fragment contains only part of the modification range: its mass is unknown
        }
      }
      if (split && span_fragments == nullptr) continue;

      Fragment f;
      f.ion_type = t;
      f.length = len;
      f.bond = prefix ? len : length - len;
      f.mono_mass = prefix_mass[last + 1] - prefix_mass[first] + ion_shifts_[t] + mod_mass;
      if (f.mono_mass <= 0 || f.mono_mass > max_mass_) continue;
      f.max_charge = (int)std::floor(f.mono_mass / min_mass_per_charge_);
      if (prsm.precursor_charge > 0) f.max_charge = std::min(f.max_charge, prsm.precursor_charge);
      if (f.max_charge < 1) continue;
      (split ? *span_fragments : fragments).push_back(f);
    }
  }
  return true;
}

bool FLASHAnnotatorAlgorithm::matchEnvelope_(const PeakArrays& spec, const PeakOccupancy& occupancy, double mono_mass, int bin, int z, double cal_factor,
                                             EnvelopeMatch& match) const
{
  const int begin = table_.offset[bin];
  if (requiredPeaks_(table_.offset[bin + 1] - begin) < 0) return false;
  // predicted m/z on the observed scale (observed / cal_factor is the calibrated m/z), and the tolerance with it
  const double scale = cal_factor / z;
  const double base = mono_mass + z * Constants::PROTON_MASS_U;
  const double tol = tol_ * 1e-6;
  // the most abundant peak is required: not searched if its window holds no peak (most envelopes)
  const double target = (base + table_.iso_index[begin + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U) * scale;
  if (! occupancy.mayHavePeak(target - target * tol, target + target * tol)) return false;
  return searchEnvelope_(spec, occupancy, base, scale, bin, match);
}

bool FLASHAnnotatorAlgorithm::searchEnvelope_(const PeakArrays& spec, const PeakOccupancy& occupancy, double base, double scale, int bin,
                                              EnvelopeMatch& match) const
{
  const int begin = table_.offset[bin];
  const int end = table_.offset[bin + 1];
  const int predicted = end - begin;
  const int required = requiredPeaks_(predicted);
  const int size = (int)spec.mzs.size();
  const int apex = begin + table_.apex_pos[bin];
  const double tol = tol_ * 1e-6;

  // the most abundant peak first: it is required
  double target = (base + table_.iso_index[apex] * Constants::C13C12_MASSDIFF_U) * scale;
  int apex_lb;
  const int apex_peak = nearestPeak(spec, 0, size, target, target * tol, apex_lb);
  if (apex_peak < 0) return false;
  const double apex_low = target - target * tol;

  match.peaks.clear();
  match.peaks.push_back(apex_peak);
  // intensities relative to the observed most abundant peak (cosine and Pearson are scale free)
  const double apex_intensity = spec.intensities[apex_peak];
  double sxy = table_.iso_rel[apex], sy = 1.0, syy = 1.0, intensity = apex_intensity;
  int matched = 1;
  int unsearched = predicted - 1;

  // isotope peaks above the most abundant one, in increasing m/z. A window that holds no peak is not searched (the search finds
  // none), so lb, where the next search starts, does not move up to it: caught up if a later window starts below a skipped one
  int lb = apex_lb;
  bool skipped = false;
  double skipped_low = 0; // the highest start of the windows skipped since the last search
  for (int j = apex + 1; j < end; j++)
  {
    target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) * scale;
    const double w = target * tol;
    int peak = -1;
    if (occupancy.mayHavePeak(target - w, target + w))
    {
      if (skipped && target - w < skipped_low) lb = std::max(lb, firstPeakAtOrAbove(spec, skipped_low));
      skipped = false;
      peak = nearestPeak(spec, lb, size, target, w, lb);
    }
    else if (! skipped || target - w > skipped_low)
    {
      skipped = true;
      skipped_low = target - w;
    }
    unsearched--;
    if (peak < 0)
    {
      if (matched + unsearched < required) return false;
      continue;
    }
    const double y = spec.intensities[peak] / apex_intensity;
    match.peaks.push_back(peak);
    matched++;
    sxy += table_.iso_rel[j] * y;
    sy += y;
    syy += y * y;
    intensity += spec.intensities[peak];
  }
  // and below it, in decreasing m/z. A window that holds no peak is skipped only if no earlier window started below it (its search
  // would start at its first peak and find none); lb, where the next search starts, is caught up as above
  lb = apex_lb;
  skipped = false; // skipped_low: the lowest start of the windows skipped since the last search
  double min_low = apex_low; // the lowest start of the windows so far
  for (int j = apex - 1; j >= begin; j--)
  {
    target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) * scale;
    const double w = target * tol;
    const double low = target - w;
    int peak = -1;
    if (low <= min_low && ! occupancy.mayHavePeak(low, target + w))
    {
      if (! skipped || low < skipped_low) skipped_low = low;
      skipped = true;
    }
    else
    {
      if (skipped && low > skipped_low) lb = std::min(lb, firstPeakAtOrAbove(spec, skipped_low));
      skipped = false;
      peak = nearestPeak(spec, 0, lb, target, w, lb);
    }
    min_low = std::min(min_low, low);
    unsearched--;
    if (peak < 0)
    {
      if (matched + unsearched < required) return false;
      continue;
    }
    const double y = spec.intensities[peak] / apex_intensity;
    match.peaks.push_back(peak);
    matched++;
    sxy += table_.iso_rel[j] * y;
    sy += y;
    syy += y * y;
    intensity += spec.intensities[peak];
  }
  if (matched < required) return false;

  const double sx = table_.rel_sum[bin];
  const double sxx = table_.rel_sum_sq[bin];
  const double cosine = sxy / std::sqrt(sxx * syy);
  if (cosine < min_cos_) return false;

  const double n = predicted;
  const double var = (n * sxx - sx * sx) * (n * syy - sy * sy);
  const double pearson = var > 0 ? (n * sxy - sx * sy) / std::sqrt(var) : 0.0;
  if (predicted > 2 && pearson < min_pearson_) return false; // two points always correlate perfectly

  match.matched = matched;
  match.apex_peak = apex_peak;
  match.cosine = cosine;
  match.pearson = pearson;
  match.intensity = intensity;
  return true;
}

double FLASHAnnotatorAlgorithm::calibrate_(const PeakArrays& spec, const std::vector<Fragment>& fragments, int& count) const
{
  count = 0;
  if (spec.mzs.empty()) return 0;
  const int size = (int)spec.mzs.size();
  const double min_mz = spec.mzs.front();
  const double max_mz = spec.mzs.back();
  const double tol = cal_tol_ * 1e-6;
  std::vector<double> errors;

  for (const auto& f : fragments)
  {
    if (! ion_for_calibration_[f.ion_type]) continue;
    const int bin = massToBin_(f.mono_mass);
    const int begin = table_.offset[bin];
    const int end = table_.offset[bin + 1];
    const int predicted = end - begin;
    if (predicted < min_peaks_) continue;
    const int apex = begin + table_.apex_pos[bin];
    const int mono = table_.mono_pos[bin] < 0 ? -1 : begin + table_.mono_pos[bin];

    for (int z = 1; z <= f.max_charge; z++)
    {
      const double base = f.mono_mass + z * Constants::PROTON_MASS_U;
      const double low = (base + table_.iso_index[begin] * Constants::C13C12_MASSDIFF_U) / z;
      const double high = (base + table_.iso_index[end - 1] * Constants::C13C12_MASSDIFF_U) / z;
      if (low * (1 - tol) > max_mz) continue;
      if (high * (1 + tol) < min_mz) break; // higher charges are further below

      // the error is measured on the monoisotopic peak if it is kept and matched, else on the most abundant one
      double target = (base + table_.iso_index[apex] * Constants::C13C12_MASSDIFF_U) / z;
      int apex_lb, lb;
      const int apex_peak = nearestPeak(spec, 0, size, target, target * tol, apex_lb);
      int mono_peak = -1;
      if (mono >= 0 && mono != apex)
      {
        const double mono_target = base / z;
        mono_peak = nearestPeak(spec, 0, apex_lb, mono_target, mono_target * tol, lb);
      }
      else if (mono == apex) { mono_peak = apex_peak; }
      if (apex_peak < 0 && mono_peak < 0) continue;

      int matched = (apex_peak >= 0 ? 1 : 0) + (mono >= 0 && mono != apex && mono_peak >= 0 ? 1 : 0);
      int unsearched = predicted - 1 - (mono >= 0 && mono != apex ? 1 : 0);
      bool fail = matched + unsearched < min_peaks_;
      lb = apex_lb;
      for (int j = apex + 1; j < end && ! fail; j++)
      {
        target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) / z;
        unsearched--;
        if (nearestPeak(spec, lb, size, target, target * tol, lb) >= 0) matched++;
        else fail = matched + unsearched < min_peaks_;
      }
      lb = apex_lb;
      for (int j = apex - 1; j >= begin && ! fail; j--)
      {
        if (j == mono) continue;
        target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) / z;
        unsearched--;
        if (nearestPeak(spec, 0, lb, target, target * tol, lb) >= 0) matched++;
        else fail = matched + unsearched < min_peaks_;
      }
      if (fail || matched < min_peaks_) continue;

      const int ref = mono_peak >= 0 ? mono : apex;
      const int ref_peak = mono_peak >= 0 ? mono_peak : apex_peak;
      const double expected = (base + table_.iso_index[ref] * Constants::C13C12_MASSDIFF_U) / z;
      errors.push_back((spec.mzs[ref_peak] - expected) / expected * 1e6);
    }
  }

  count = (int)errors.size();
  if (count < min_calibration_envelopes_) return 0;
  const Size mid = errors.size() / 2;
  std::nth_element(errors.begin(), errors.begin() + mid, errors.end());
  double median = errors[mid];
  if (errors.size() % 2 == 0)
  {
    median = (median + *std::max_element(errors.begin(), errors.begin() + mid)) / 2.0;
  }
  return median;
}

void FLASHAnnotatorAlgorithm::annotateOne_(const PeakArrays& spec, const std::vector<double>& sorted_intensities, const PrSM& prsm, Result& result) const
{
  auto density_at = [&spec](double mz) { return densityAt_(spec, mz); };
  auto rank_of = [&sorted_intensities](double intensity) { return rankOf_(sorted_intensities, intensity); };
  std::vector<Fragment> fragments, span_fragments;
  if (! getFragments(prsm, fragments, &span_fragments)) return;

  const int length = (int)prsm.sequence.size();
  result.annotated = true;
  result.envelopes.clear();

  int cal_count = 0;
  const double cal_ppm = calibrate_on_ ? calibrate_(spec, fragments, cal_count) : 0.0;
  const double cal_factor = 1.0 + cal_ppm * 1e-6;
  result.calibration_ppm = cal_ppm;
  result.calibration_envelopes = cal_count;

  const bool empty = spec.mzs.empty();
  const double min_mz = empty ? 0 : spec.mzs.front();
  const double max_mz = empty ? 0 : spec.mzs.back();
  const double tol = tol_ * 1e-6;
  PeakOccupancy occupancy;
  occupancy.build(spec);

  // fragment hypotheses (ion type group, length): ion types whose envelopes fall on the same peaks count once
  const int hypothesis_count = ion_group_count_ * std::max(1, length);
  std::vector<char> testable(hypothesis_count, 0), matched(hypothesis_count, 0);

  // target fragments: every accepted envelope is recorded
  std::vector<char> covered(std::max(1, length), 0);
  std::vector<char> terminal_peaks(spec.mzs.size(), 0); // peaks of accepted terminal envelopes
  EnvelopeMatch match;
  for (const auto& f : fragments)
  {
    if (empty) break; // nothing can match; no hypothesis is testable
    const int bin = massToBin_(f.mono_mass);
    const int predicted = table_.offset[bin + 1] - table_.offset[bin];
    if (requiredPeaks_(predicted) < 0) continue; // can never be accepted, not a hypothesis
    const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
    const int hypothesis = ion_group_[f.ion_type] * length + f.length;
    for (int z = 1; z <= f.max_charge; z++)
    {
      const double apex_mz = (f.mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z;
      const double observed = apex_mz * cal_factor;
      if (observed * (1 - tol) > max_mz) continue;
      if (observed * (1 + tol) < min_mz) break; // higher charges are further below
      testable[hypothesis] = 1;
      if (! matchEnvelope_(spec, occupancy, f.mono_mass, bin, z, cal_factor, match)) continue;
      matched[hypothesis] = 1;
      covered[f.bond] = 1;
      for (int k : match.peaks) terminal_peaks[k] = 1;

      Envelope e;
      e.ion_type = f.ion_type;
      e.length = f.length;
      e.bond = f.bond;
      e.charge = z;
      e.mono_mass = f.mono_mass;
      e.most_abundant_mz = apex_mz;
      e.matched_peaks = match.matched;
      e.predicted_peaks = predicted;
      e.cosine = match.cosine;
      e.pearson = match.pearson;
      e.intensity = match.intensity;
      e.ppm_error = (spec.mzs[match.apex_peak] / cal_factor - apex_mz) / apex_mz * 1e6;
      e.intensity_rank = rank_of(spec.intensities[match.apex_peak]);
      e.local_density = density_at(spec.mzs[match.apex_peak]);
      result.envelopes.push_back(e);
    }
  }
  // fragments across an unlocalized modification range: their mass is unknown, so they are never matched, but they are hypotheses the
  // proteoform makes (a long range would otherwise shrink the hypotheses to the few fragments outside it and inflate the score)
  for (const auto& f : span_fragments)
  {
    if (empty) break;
    const int bin = massToBin_(f.mono_mass);
    if (requiredPeaks_(table_.offset[bin + 1] - table_.offset[bin]) < 0) continue;
    const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
    for (int z = 1; z <= f.max_charge; z++)
    {
      const double observed = (f.mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z * cal_factor;
      if (observed * (1 - tol) > max_mz) continue;
      if (observed * (1 + tol) < min_mz) break;
      testable[ion_group_[f.ion_type] * length + f.length] = 1;
      break;
    }
  }
  const int hypotheses = (int)std::count(testable.begin(), testable.end(), 1);
  const int matched_hypotheses = (int)std::count(matched.begin(), matched.end(), 1);

  // decoy fragments: each target fragment (its isotope pattern and charges) moved by a mass shift; counts only, and only
  // for the testable target hypotheses, so that at most K * n decoy hypotheses match
  int decoy_matched = 0;
  std::vector<char> decoy_matched_hypothesis(hypothesis_count, 0);
  for (Size s = 0; s < decoy_shifts_.size() && hypotheses > 0; s++)
  {
    const double shift = decoy_shifts_[s];
    std::fill(decoy_matched_hypothesis.begin(), decoy_matched_hypothesis.end(), 0);
    for (Size i = 0; i < fragments.size() + span_fragments.size(); i++)
    {
      const auto& f = i < fragments.size() ? fragments[i] : span_fragments[i - fragments.size()];
      const int hypothesis = ion_group_[f.ion_type] * length + f.length;
      if (! testable[hypothesis] || (decoy_matched_hypothesis[hypothesis] && ! report_decoys_)) continue;
      const double mono_mass = f.mono_mass + shift;
      if (mono_mass <= 0) continue;
      const int bin = massToBin_(f.mono_mass);
      if (requiredPeaks_(table_.offset[bin + 1] - table_.offset[bin]) < 0) continue;
      const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
      for (int z = 1; z <= f.max_charge; z++)
      {
        const double apex_mz = (mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z;
        const double observed = apex_mz * cal_factor;
        if (observed * (1 - tol) > max_mz) continue;
        if (observed * (1 + tol) < min_mz) break;
        if (matchEnvelope_(spec, occupancy, mono_mass, bin, z, cal_factor, match))
        {
          if (! decoy_matched_hypothesis[hypothesis]) decoy_matched++; // a hypothesis counts once, whatever the number of charges and ion types
          decoy_matched_hypothesis[hypothesis] = 1;
          if (! report_decoys_) break;
          Envelope e;
          e.ion_type = f.ion_type;
          e.length = f.length;
          e.bond = f.bond;
          e.charge = z;
          e.mono_mass = mono_mass;
          e.most_abundant_mz = apex_mz;
          e.matched_peaks = match.matched;
          e.predicted_peaks = table_.offset[bin + 1] - table_.offset[bin];
          e.cosine = match.cosine;
          e.pearson = match.pearson;
          e.intensity = match.intensity;
          e.ppm_error = (spec.mzs[match.apex_peak] / cal_factor - apex_mz) / apex_mz * 1e6;
          e.intensity_rank = rank_of(spec.intensities[match.apex_peak]);
          e.local_density = density_at(spec.mzs[match.apex_peak]);
          e.decoy_shift = shift;
          result.envelopes.push_back(e);
        }
      }
    }
  }

  int covered_bonds = 0;
  for (int b = 1; b < length; b++)
  {
    covered_bonds += covered[b];
  }
  result.covered_bonds = covered_bonds;
  result.coverage = length > 1 ? 100.0 * covered_bonds / (length - 1) : 0.0;
  result.matched_fragments = matched_hypotheses;
  result.fragment_hypotheses = hypotheses;
  result.decoy_matched_fragments = decoy_matched;

  if (! internal_on_ || empty || length < 5) return;

  // internal pieces [first, last], both ends inside the proteoform, 2 .. internal_max_length residues, every internal ion type
  const auto& aa = residueMasses();
  std::vector<double> prefix(length + 1, .0);
  for (int i = 0; i < length; i++)
  {
    prefix[i + 1] = prefix[i] + aa[(unsigned char)prsm.sequence[i]];
  }
  InternalPieces pieces;
  pieces.length = length;
  pieces.offsets = (int)internal_offsets_.size();
  pieces.max_length = internal_max_length_;
  pieces.base.assign(length, 0);
  int total = 0;
  for (int first = 1; first < length; first++)
  {
    pieces.base[first] = total;
    if (first <= length - 3) total += (pieces.lastMax(first) - first) * pieces.offsets;
  }
  pieces.mono_mass.assign(total, -1.0);
  pieces.max_charge.assign(total, 0);
  int testable_pieces = 0; // pieces (first, last) testable with at least one internal ion type
  for (int first = 1; first <= length - 3; first++)
  {
    for (int last = first + 1; last <= pieces.lastMax(first); last++)
    {
      double mass = prefix[last + 1] - prefix[first];
      bool split = false;
      for (const auto& mod : prsm.modifications)
      {
        const int ms = std::min(mod.start, mod.end), me = std::max(mod.start, mod.end);
        if (ms >= first && me <= last) mass += mod.mass;
        else if (me >= first && ms <= last) split = true; // the piece splits the modification range: its mass is unknown
      }
      if (split) continue;
      bool testable = false;
      for (int o = 0; o < pieces.offsets; o++)
      {
        const double mono_mass = mass + internal_offsets_[o];
        if (mono_mass <= 0 || mono_mass > max_mass_) continue;
        int max_charge = (int)std::floor(mono_mass / min_mass_per_charge_);
        if (prsm.precursor_charge > 0) max_charge = std::min(max_charge, prsm.precursor_charge);
        if (max_charge < 1) continue;
        const int bin = massToBin_(mono_mass);
        if (table_.offset[bin + 1] - table_.offset[bin] < min_peaks_) continue;
        // testable: a charge whose most abundant peak lies in the spectrum's m/z range
        const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
        bool in_range = false;
        for (int z = 1; z <= max_charge && ! in_range; z++)
        {
          const double observed = (mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z * cal_factor;
          in_range = observed * (1 - tol) <= max_mz && observed * (1 + tol) >= min_mz;
        }
        if (! in_range) continue;
        const int idx = pieces.index(first, last, o);
        pieces.mono_mass[idx] = mono_mass;
        pieces.max_charge[idx] = max_charge;
        testable = true;
      }
      if (testable) testable_pieces++;
    }
  }

  // no extension where the proteoform's start or end is undetermined: its "internal" pieces may touch the real terminus
  const int extension_min_run = (prsm.known_n_term && prsm.known_c_term) ? extension_min_run_ : 0;
  std::vector<char> authentic, decoy_authentic;
  const auto counted = annotateInternal_(spec, occupancy, sorted_intensities, pieces, prsm.sequence, terminal_peaks, cal_factor, 0.0, true,
                                         extension_min_run, std::vector<char>(length, 0), result.envelopes, authentic);
  // the decoy worlds of the internal fragments are only needed to study chance (the internal p-value comes from the calibration
  // model): they run with report_decoy_envelopes only
  int decoy_counted = -1;
  for (Size s = 0; report_decoys_ && s < decoy_shifts_.size(); s++)
  {
    decoy_counted = std::max(0, decoy_counted);
    decoy_counted += (int)annotateInternal_(spec, occupancy, sorted_intensities, pieces, prsm.sequence, terminal_peaks, cal_factor, decoy_shifts_[s],
                                            true, extension_min_run, authentic, result.envelopes, decoy_authentic)
                       .size();
  }

  std::vector<char> covered_all = covered;
  for (const auto& [first, last] : counted)
  {
    covered_all[first] = 1;    // bond between first - 1 and first
    covered_all[last + 1] = 1; // bond between last and last + 1
  }
  int covered_all_bonds = 0;
  for (int b = 1; b < length; b++)
  {
    covered_all_bonds += covered_all[b];
  }
  result.internal_hypotheses = testable_pieces;
  result.internal_matched = (int)counted.size();
  result.internal_decoy_matched = decoy_counted;
  result.authentic_cleavages = (int)std::count(authentic.begin(), authentic.end(), 1);
  result.covered_bonds_with_internal = covered_all_bonds;
  result.coverage_with_internal = 100.0 * covered_all_bonds / (length - 1);
}

std::set<std::pair<int, int>> FLASHAnnotatorAlgorithm::annotateInternal_(const PeakArrays& spec, const PeakOccupancy& occupancy,
                                                                         const std::vector<double>& sorted_intensities,
                                                                         const InternalPieces& pieces, const std::string& sequence,
                                                                         const std::vector<char>& terminal_peaks, double cal_factor, double shift,
                                                                         bool record, int extension_min_run, const std::vector<char>& given_authentic,
                                                                         std::vector<Envelope>& envelopes, std::vector<char>& authentic) const
{
  const double tol = tol_ * 1e-6;
  const double min_mz = spec.mzs.front();
  const double max_mz = spec.mzs.back();
  const int length = pieces.length;
  const int offsets = pieces.offsets;

  // every piece at every charge; the envelopes of accepted pieces are kept for the record
  std::vector<char> accepted(pieces.mono_mass.size(), 0);
  std::unordered_map<int, std::vector<Envelope>> piece_envelopes;
  EnvelopeMatch match;
  for (int first = 1; first <= length - 3; first++)
  {
    for (int last = first + 1; last <= pieces.lastMax(first); last++)
    {
      for (int o = 0; o < offsets; o++)
      {
        const int idx = pieces.index(first, last, o);
        if (pieces.mono_mass[idx] < 0) continue;
        const double mono_mass = pieces.mono_mass[idx] + shift;
        if (mono_mass <= 0) continue;
        const int bin = massToBin_(pieces.mono_mass[idx]);
        const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
        for (int z = 1; z <= pieces.max_charge[idx]; z++)
        {
          const double apex_mz = (mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z;
          const double observed = apex_mz * cal_factor;
          if (observed * (1 - tol) > max_mz) continue;
          if (observed * (1 + tol) < min_mz) break;
          if (! matchEnvelope_(spec, occupancy, mono_mass, bin, z, cal_factor, match) || terminal_peaks[match.apex_peak]) continue;
          accepted[idx] = 1;
          if (! record) break;
          Envelope e;
          e.ion_type = o;
          e.first = first;
          e.length = last - first + 1;
          e.bond = first;
          e.charge = z;
          e.mono_mass = mono_mass;
          e.most_abundant_mz = apex_mz;
          e.matched_peaks = match.matched;
          e.predicted_peaks = table_.offset[bin + 1] - table_.offset[bin];
          e.cosine = match.cosine;
          e.pearson = match.pearson;
          e.intensity = match.intensity;
          e.ppm_error = (spec.mzs[match.apex_peak] / cal_factor - apex_mz) / apex_mz * 1e6;
          e.intensity_rank = rankOf_(sorted_intensities, spec.intensities[match.apex_peak]);
          e.local_density = densityAt_(spec, spec.mzs[match.apex_peak]);
          e.decoy_shift = shift;
          piece_envelopes[idx].push_back(e);
        }
      }
    }
  }

  // ladders: for the cleavage after residue a and one internal ion type, accepted pieces of consecutive lengths that all end with
  // a, or all start right after it. run_end / run_start: the size of a piece's run in the group of its last / first residue.
  authentic.assign(length, 0);
  // threshold by the residues around the cleavage after residue a: after Asp, before Pro (the smaller one where both hold), other
  auto min_run = [&sequence, length, this](int a) {
    const bool asp = sequence[a] == 'D' || sequence[a] == 'd';
    const bool pro = a + 1 < length && (sequence[a + 1] == 'P' || sequence[a + 1] == 'p');
    if (asp && pro) return std::min(min_run_asp_, min_run_pro_);
    if (asp) return min_run_asp_;
    return pro ? min_run_pro_ : min_run_other_;
  };
  std::vector<int> run_end(accepted.size(), 0), run_start(accepted.size(), 0);
  std::vector<int> run;
  auto flush = [&](int a, std::vector<int>& run_size) {
    if ((int)run.size() >= min_run(a)) authentic[a] = 1;
    for (int k : run) run_size[k] = (int)run.size();
    run.clear();
  };
  for (int a = 0; a <= length - 2; a++)
  {
    for (int o = 0; o < offsets; o++)
    {
      // pieces [first, a] in order of increasing length
      int prev_first = -1;
      for (int first = a - 1; a <= length - 2 && first >= 1 && a - first < internal_max_length_; first--)
      {
        const int idx = pieces.index(first, a, o);
        if (! accepted[idx] || prev_first != first + 1) flush(a, run_end);
        if (accepted[idx])
        {
          run.push_back(idx);
          prev_first = first;
        }
        else prev_first = -1;
      }
      flush(a, run_end);
      // pieces [a + 1, last] in order of increasing length
      int prev_last = -1;
      for (int last = a + 2; a + 1 <= length - 3 && last <= pieces.lastMax(a + 1); last++)
      {
        const int idx = pieces.index(a + 1, last, o);
        if (! accepted[idx] || prev_last != last - 1) flush(a, run_start);
        if (accepted[idx])
        {
          run.push_back(idx);
          prev_last = last;
        }
        else prev_last = -1;
      }
      flush(a, run_start);
    }
  }

  std::set<std::pair<int, int>> counted;
  for (int first = 1; first <= length - 3; first++)
  {
    for (int last = first + 1; last <= pieces.lastMax(first); last++)
    {
      for (int o = 0; o < offsets; o++)
      {
        const int idx = pieces.index(first, last, o);
        if (! accepted[idx]) continue;
        const int a_end = last, a_start = first - 1;
        bool count = run_end[idx] >= min_run(a_end) || run_start[idx] >= min_run(a_start); // ladder pieces
        if (! count && extension_min_run > 0)
        {
          count = ((authentic[a_end] || given_authentic[a_end]) && run_end[idx] >= extension_min_run)
                  || ((authentic[a_start] || given_authentic[a_start]) && run_start[idx] >= extension_min_run);
        }
        if (! count) continue;
        counted.emplace(first, last);
        if (record)
        {
          const auto it = piece_envelopes.find(idx);
          if (it != piece_envelopes.end()) envelopes.insert(envelopes.end(), it->second.begin(), it->second.end());
        }
      }
    }
  }
  return counted;
}

double FLASHAnnotatorAlgorithm::binomialTailScore(int n, int k, double p)
{
  if (k <= 0 || n <= 0) return 0;
  if (k > n) k = n;
  if (p <= 0) return max_score_;
  if (p >= 1) return 0;

  const double log_p = std::log(p);
  const double log_q = std::log1p(-p);
  const double log_n_fact = std::lgamma(n + 1.0);
  auto log_pmf = [&](int j) { return log_n_fact - std::lgamma(j + 1.0) - std::lgamma(n - j + 1.0) + j * log_p + (n - j) * log_q; };

  // the terms of the upper tail grow up to the mode and decrease after it; sum them relative to the largest one
  const int mode = std::min(n, (int)std::floor((n + 1) * p));
  const int top = std::max(k, mode);
  const double log_top = log_pmf(top);
  double sum = 0;
  for (int j = k; j <= n; j++)
  {
    const double d = log_pmf(j) - log_top;
    if (j > top && d < -50) break; // the rest is negligible
    sum += std::exp(d);
  }
  const double score = -(log_top + std::log(sum)) / std::log(10.0);
  return std::max(0.0, std::min(max_score_, score));
}

void FLASHAnnotatorAlgorithm::annotate(const std::vector<PeakArrays>& spectra, const std::vector<PrSM>& prsms, std::vector<Result>& results)
{
  results.assign(prsms.size(), Result());

  // the isotope table must cover the heaviest fragment. It is computed here, once, before the parallel region.
  const auto& aa = residueMasses();
  double max_fragment_mass = 0;
  for (const auto& prsm : prsms)
  {
    double mass = 0;
    for (const char c : prsm.sequence)
    {
      mass += aa[(unsigned char)c];
    }
    for (const auto& mod : prsm.modifications)
    {
      mass += std::max(.0, mod.mass);
    }
    max_fragment_mass = std::max(max_fragment_mass, mass + 100.0); // + the largest ion type shift
  }
  ensureIsotopeTable_(max_fragment_mass);

  // longest proteoforms first, so that the dynamic schedule does not end with a long one
  std::vector<int> order(prsms.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&prsms](int a, int b) { return prsms[a].sequence.size() > prsms[b].sequence.size(); });

  const int spectrum_count = (int)spectra.size();
  // sorted intensities of the spectra that are annotated (for the intensity rank of the envelopes)
  std::vector<std::vector<double>> sorted_intensities(spectra.size());
  std::vector<char> used(spectra.size(), 0);
  for (const auto& prsm : prsms)
  {
    if (prsm.spectrum_index >= 0 && prsm.spectrum_index < spectrum_count) used[prsm.spectrum_index] = 1;
  }
#pragma omp parallel for schedule(dynamic)
  for (int s = 0; s < spectrum_count; s++)
  {
    if (! used[s]) continue;
    sorted_intensities[s] = spectra[s].intensities;
    std::sort(sorted_intensities[s].begin(), sorted_intensities[s].end());
  }

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int)order.size(); i++)
  {
    const int index = order[i];
    const PrSM& prsm = prsms[index];
    if (prsm.spectrum_index < 0 || prsm.spectrum_index >= spectrum_count) continue;
    annotateOne_(spectra[prsm.spectrum_index], sorted_intensities[prsm.spectrum_index], prsm, results[index]);
  }

  // scores (serial: std::lgamma need not be thread safe). Without a testable hypothesis there is nothing to score.
  // The chance rate is that of the spectrum: the decoy fragments of every PrSM annotated on it here, over all their hypotheses. A rate of
  // each PrSM's own decoys (a few matches among some hundred tests) differed several-fold between PrSMs of one spectrum by chance, which
  // decided between them more than their matched fragments did.
  const double decoy_sets = (double)decoy_shifts_.size();
  std::vector<double> spectrum_decoys(spectrum_count, 0.0), spectrum_tests(spectrum_count, 0.0);
  for (Size i = 0; i < results.size(); i++)
  {
    const auto& r = results[i];
    if (! r.annotated || decoy_shifts_.empty() || r.fragment_hypotheses <= 0) continue;
    spectrum_decoys[prsms[i].spectrum_index] += r.decoy_matched_fragments;
    spectrum_tests[prsms[i].spectrum_index] += decoy_sets * r.fragment_hypotheses;
  }
  for (Size i = 0; i < results.size(); i++)
  {
    auto& r = results[i];
    if (! r.annotated || decoy_shifts_.empty() || r.fragment_hypotheses <= 0) continue;
    const int n = r.fragment_hypotheses;
    const int k = r.matched_fragments;
    const int s = prsms[i].spectrum_index;
    const double p0 = (spectrum_decoys[s] + 1.0) / (spectrum_tests[s] + 2.0);
    r.chance_rate = p0;
    r.score = binomialTailScore(n, k, p0);
    r.llr = k * std::log(p1_ / p0) + (n - k) * std::log((1.0 - p1_) / (1.0 - p0));
    if (r.internal_hypotheses > 0)
    {
      // from the calibration model: the chance of at least one counted piece, and -log10 P(K >= k | n)
      r.internal_chance_rate = std::pow(10.0, -internalLog10P_(1, r.internal_hypotheses));
      r.internal_score = internalLog10P_(r.internal_matched, r.internal_hypotheses);
    }
    // FAScore: the calibrated terminal and internal p-values combined (Fisher, chi-square with 4 degrees of freedom)
    const double log10_pt = calibratedTerminalLog10P_(r.score);
    const double log10_pi = internalLog10P_(r.internal_matched, r.internal_hypotheses);
    const double x = 2.0 * std::log(10.0) * (log10_pt + log10_pi);           // -2 ln(p_t p_i)
    const double log_p = -x / 2.0 + std::log1p(x / 2.0);                     // ln P(chi2_4 >= x)
    r.fa_score = std::min(max_score_, std::max(0.0, -log_p / std::log(10.0)));
  }
}

namespace
{
  /// keeps, of the spectra read by MzMLFile::transform, the first MSn spectrum of every needed scan number
  class NeededSpectraConsumer : public Interfaces::IMSDataConsumer
  {
  public:
    /// with @p scan_from_index, scan numbers are spectrum index + 1 unless the native ID carries one (all spectra must be read)
    NeededSpectraConsumer(const std::set<int>& needed_scans, bool scan_from_index, std::map<int, MSSpectrum>& out):
      needed_scans_(needed_scans), scan_from_index_(scan_from_index), out_(out)
    {
    }
    void setExperimentalSettings(const ExperimentalSettings& settings) override
    {
      // as FLASHDeconvAlgorithm::getScanNumber
      type_accession_ = "MS:1000768";
      if (! settings.getSourceFiles().empty() && ! settings.getSourceFiles()[0].getNativeIDTypeAccession().empty())
      {
        type_accession_ = settings.getSourceFiles()[0].getNativeIDTypeAccession();
      }
    }
    void setExpectedSize(Size, Size) override
    {
    }
    void consumeSpectrum(SpectrumType& spec) override
    {
      const int index = read_count_++;
      std::vector<std::string> native_ids;
      StringUtils::split(spec.getNativeID(), ",", native_ids);
      int scan = native_ids.empty() ? -1 : SpectrumLookup::extractScanNumber(native_ids.back(), type_accession_);
      if (scan < 0)
      {
        without_scan_number_ = true;
        if (! scan_from_index_) return;
        scan = index + 1;
      }
      if (spec.getMSLevel() < 2 || needed_scans_.find(scan) == needed_scans_.end() || out_.find(scan) != out_.end()) return;
      out_[scan] = std::move(spec);
    }
    void consumeChromatogram(ChromatogramType&) override
    {
    }
    int readCount() const
    {
      return read_count_;
    }
    bool withoutScanNumber() const
    {
      return without_scan_number_;
    }

  private:
    const std::set<int>& needed_scans_;
    bool scan_from_index_;
    std::map<int, MSSpectrum>& out_;
    std::string type_accession_ = "MS:1000768";
    int read_count_ = 0;
    bool without_scan_number_ = false;
  };
} // namespace

Size FLASHAnnotatorAlgorithm::loadSpectra(const std::string& mzml_file, const std::set<int>& scans, std::map<int, MSSpectrum>& spectra)
{
  if (FileHandler::getTypeByFileName(mzml_file) == FileTypes::RAW)
  {
#ifdef WITH_THERMO_RAW
    // a Thermo .raw file, read with Thermo's peak picking (centroided) as FLASHHelperClasses::loadExperiment reads it, MSn spectra only;
    // its native IDs always hold the scan number
    PeakMap map;
    FileHandler file_handler;
    PeakFileOptions options = file_handler.getOptions();
    options.setMSLevels({2, 3, 4, 5});
    file_handler.setOptions(options);
    file_handler.loadExperiment(mzml_file, map, {FileTypes::RAW}, ProgressLogger::NONE);
    NeededSpectraConsumer consumer(scans, false, spectra);
    for (auto& spec : map) consumer.consumeSpectrum(spec);
    return consumer.readCount();
#else
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "'" + mzml_file + "' is a Thermo .raw file, which this build cannot read (OpenMS built without WITH_THERMO_RAW): "
                                      "convert it to mzML first, e.g. with ThermoRawFileParser.");
#endif
  }
  MzMLFile mzml;
  mzml.setLogType(ProgressLogger::NONE);
  PeakFileOptions options = mzml.getOptions();
  options.setMSLevels({2, 3, 4, 5});
  // chromatograms are not skipped: with setSkipChromatograms(true), MzMLHandler ends a progress it did not start
  mzml.setOptions(options);
  NeededSpectraConsumer consumer(scans, false, spectra);
  mzml.transform(mzml_file, &consumer, true);
  Size read = consumer.readCount();
  if (consumer.withoutScanNumber())
  {
    // scan numbers then come from the spectrum indices, which require all spectra
    spectra.clear();
    options.setMSLevels({});
    mzml.setOptions(options);
    NeededSpectraConsumer all(scans, true, spectra);
    mzml.transform(mzml_file, &all, true);
    read = all.readCount();
  }
  return read;
}

// At most this many peaks per spectrum (FLASHDeconv's cap for centroided spectra): a spectrum with more (e.g. a long-transient profile
// spectrum, 327,756 centroided peaks) keeps its most intense ones. In such a spectrum about 17 % of random fragment hypotheses found an
// envelope (0.4-2 % in the others), shuffled sequences reached 32 % terminal and 100 % internal coverage; with the cap 8 % and 12 %.
static constexpr Size max_peak_count_ = 30000;

static void keepMostIntensePeaks_(FLASHAnnotatorAlgorithm::PeakArrays& peaks)
{
  const Size n = peaks.mzs.size();
  if (n <= max_peak_count_) return;
  std::vector<Size> order(n);
  std::iota(order.begin(), order.end(), 0);
  // the most intense first, equal intensities by m/z (deterministic)
  std::nth_element(order.begin(), order.begin() + max_peak_count_, order.end(), [&](Size l, Size r) {
    if (peaks.intensities[l] != peaks.intensities[r]) return peaks.intensities[l] > peaks.intensities[r];
    return l < r;
  });
  std::vector<char> keep(n, 0);
  for (Size k = 0; k < max_peak_count_; k++) keep[order[k]] = 1;
  Size o = 0;
  for (Size t = 0; t < n; t++)
  {
    if (! keep[t]) continue;
    peaks.mzs[o] = peaks.mzs[t];
    peaks.intensities[o] = peaks.intensities[t];
    o++;
  }
  peaks.mzs.resize(o);
  peaks.intensities.resize(o);
  peaks.indexBuckets();
}

bool FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(const MSSpectrum& spec, PeakArrays& peaks)
{
  if (spec.getType(true) == SpectrumSettings::SpectrumType::PROFILE)
  {
    PeakPickerHiRes picker;
    MSSpectrum picked;
    picker.pick(spec, picked);
    toPeakArrays(picked, peaks);
    keepMostIntensePeaks_(peaks);
    return true;
  }
  toPeakArrays(spec, peaks);
  keepMostIntensePeaks_(peaks);
  return false;
}

int FLASHAnnotatorAlgorithm::precursorChargeCap(const MSSpectrum& raw, double precursor_mass, int& reason)
{
  int raw_charge = 0;
  double precursor_mz = 0, low = 0, high = 0;
  if (! raw.getPrecursors().empty())
  {
    const auto& precursor = raw.getPrecursors()[0];
    raw_charge = precursor.getCharge();
    precursor_mz = precursor.getMZ();
    low = precursor_mz - precursor.getIsolationWindowLowerOffset();
    high = precursor_mz + precursor.getIsolationWindowUpperOffset();
  }
  std::vector<double> iso_mzs, iso_intensities;
  // true if an isotope peak of mass at charge z lies in [low, high]
  auto in_window = [&](double mass, int z) {
    getIsotopeEnvelope(mass, z, iso_mzs, iso_intensities);
    return ! iso_mzs.empty() && iso_mzs.front() <= high && iso_mzs.back() >= low;
  };
  // the charge at which the precursor mass has the isolation target m/z
  const int implied_charge = (precursor_mass > 0 && precursor_mz > Constants::PROTON_MASS_U)
                               ? std::max(1, (int)std::round(precursor_mass / (precursor_mz - Constants::PROTON_MASS_U)))
                               : 0;
  if (raw_charge > 0)
  {
    // a harmonic or co-isolated raw charge does not fit the proteoform: fragments above it would never be tried
    if (implied_charge <= 0 || implied_charge == raw_charge || in_window(precursor_mass, raw_charge))
    {
      reason = 0;
      return raw_charge;
    }
    if (in_window(precursor_mass, implied_charge))
    {
      reason = 1;
      return implied_charge;
    }
    reason = 2;
    return std::max(raw_charge, implied_charge);
  }
  if (implied_charge > 0)
  {
    reason = 3;
    return implied_charge;
  }
  reason = 4;
  return 0;
}
} // namespace OpenMS
