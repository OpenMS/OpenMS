// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include "ProSEImpl.h"

#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <functional>
#include <numeric>
#include <iomanip>
#include <locale>
#include <ostream>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  namespace
  {
    // The peaks of @p s corrected by @p model (an MzRecalibration_::AtRT), each at the m/z it was measured at: a deisotoped
    // peak of charge z > 1 (@p charges, aligned with the peaks; empty: every peak as measured) holds the m/z of its singly
    // charged ion, z mz - (z - 1) proton (Deisotoper::deisotopeAndSingleCharge()), and keeps it
    template <class Model>
    void correctAtMeasuredMZ(MSSpectrum& s, const Model& model, const std::vector<std::uint8_t>& charges)
    {
      const double proton = Constants::PROTON_MASS_U;
      for (Size k = 0; k < s.size(); ++k)
      {
        const double z = charges.empty() ? 1.0 : static_cast<double>(std::max<std::uint8_t>(1, charges[k]));
        const double mz = s[k].getMZ();
        s[k].setMZ(z == 1.0 ? model.correct(mz) : model.correct((mz + (z - 1.0) * proton) / z) * z - (z - 1.0) * proton);
      }
      if (!s.isSorted()) { s.sortByPosition(); } // (peaks of different charge may change places)
    }
  }

  // =====================================================================
  // Helper: decoy-marker membership test (prefix/suffix), pure std::string.
  // =====================================================================
  bool ProSEAlgorithm::accessionHasDecoyMarker_(const std::string& accession,
                                                const std::string& marker, bool is_prefix)
  {
    if (marker.empty() || accession.size() < marker.size()) { return false; }
    return is_prefix ? (accession.compare(0, marker.size(), marker) == 0)
                     : (accession.compare(accession.size() - marker.size(), marker.size(), marker) == 0);
  }

  bool ProSEAlgorithm::useSearchEngineScore_(PeptideIdentificationList& ids)
  {
    // The name of the score of each PSM: none for those that have it as main score already
    const std::string main_score = SEARCH_ENGINE_SCORE, meta_score = main_score + "_score";
    std::vector<const std::string*> source(ids.size(), nullptr);
    for (Size i = 0; i < ids.size(); ++i)
    {
      const PeptideIdentification& id = ids[i];
      if (id.getHits().empty() || id.getScoreType() == main_score) { continue; }
      const PeptideHit& hit = id.getHits()[0];
      if (hit.metaValueExists(main_score)) { source[i] = &main_score; }
      else if (hit.metaValueExists(meta_score)) { source[i] = &meta_score; }
      else { return false; }
    }
    // one switcher per meta value name
    const auto switcher = [&main_score](const std::string& name)
    {
      IDScoreSwitcherAlgorithm s;
      Param p = s.getParameters();
      p.setValue("new_score", name);
      p.setValue("new_score_type", main_score);
      p.setValue("new_score_orientation", "higher_better");
      s.setParameters(p);
      return s;
    };
    IDScoreSwitcherAlgorithm from_main = switcher(main_score), from_meta = switcher(meta_score);
    Size counter = 0;
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (source[i] != nullptr) { (source[i] == &main_score ? from_main : from_meta).switchScores(ids[i], counter); }
    }
    return true;
  }

  Size ProSEAlgorithm::countDecoyEntries_(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& marker, bool is_prefix)
  {
    return static_cast<Size>(std::count_if(db.begin(), db.end(), [&](const FASTAFile::FASTAEntry& e)
                                           { return accessionHasDecoyMarker_(e.identifier, marker, is_prefix); }));
  }

  // =====================================================================
  // Helper: capture pre-FDR statistics (target/decoy counts + HyperScore
  // distribution). Must run BEFORE FalseDiscoveryRate, which overwrites the
  // hit score with the q-value and may drop decoy hits.
  // =====================================================================
  void ProSEAlgorithm::capturePreFdrStats_(const PeptideIdentificationList& peptide_ids,
                                           RunStatistics& stats)
  {
    stats.target_psms = 0;
    stats.decoy_psms = 0;
    // The top hits of the spectra are read in parallel (the target/decoy value by name index: looking the name up takes
    // a lock on the registry, which the threads would contend for; by name if it is not registered), then their scores
    // collected and counted in spectrum order, as before (also where an error stops it).
    enum TopHit : char { NO_HIT, OTHER, TARGET, DECOY };
    const UInt td_index = MetaInfoInterface::metaRegistry().getIndex(Constants::UserParam::TARGET_DECOY);
    const SignedSize num_ids = static_cast<SignedSize>(peptide_ids.size());
    std::vector<TopHit> top_hits(peptide_ids.size(), NO_HIT);
    std::vector<double> top_scores(peptide_ids.size());
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(static, 1024)
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      try
      {
        const PeptideIdentification& pid = peptide_ids[i];
        if (pid.getHits().empty()) { continue; }
        const PeptideHit& top = pid.getHits().front();
        top_hits[i] = OTHER;
        top_scores[i] = top.getScore();

        const bool by_name = td_index == UInt(-1);
        if (by_name ? !top.metaValueExists(Constants::UserParam::TARGET_DECOY) : !top.metaValueExists(td_index)) { continue; }
        const std::string td = (by_name ? top.getMetaValue(Constants::UserParam::TARGET_DECOY) : top.getMetaValue(td_index)).toString();
        // "target+decoy" counts as target (OpenMS FDR semantics).
        if (td == "decoy") { top_hits[i] = DECOY; }
        else if (!td.empty()) { top_hits[i] = TARGET; }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    std::vector<double> hyperscores;
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      if (top_hits[i] == NO_HIT) { continue; }
      hyperscores.push_back(top_scores[i]);
      if (errors[i]) { std::rethrow_exception(errors[i]); }
      if (top_hits[i] == DECOY) { ++stats.decoy_psms; }
      else if (top_hits[i] == TARGET) { ++stats.target_psms; }
    }

    if (!hyperscores.empty())
    {
      auto [mn, mx] = std::minmax_element(hyperscores.begin(), hyperscores.end());
      stats.hyperscore_min = *mn;
      stats.hyperscore_max = *mx;
      stats.hyperscore_median = Math::median(hyperscores.begin(), hyperscores.end());
      stats.score_stats_valid = true;
    }
  }

  namespace
  {
    // What collectRunStatistics_() and updateFinalStats() count over the top hits
    struct TopHitSummary
    {
      Size matched = 0, n_target = 0, n_decoy = 0;
      std::map<Int, Size> charges;
      std::map<Size, Size> missed_cleavages;
      Size unique_peptides = 0, unique_proteins = 0; // distinct sequence strings and protein accessions
      std::vector<double> precursor_errors, fragment_errors; // top hits in input order (with_errors only)
      Size fragment_errors_psms = 0; // the PSMs with any fragment error (with_errors only)
    };

    // One pass over the top hits, in parallel: each thread takes a contiguous range of the identifications, and the
    // ranges are merged in order, so the error lists come out in input order as from a single pass. The regular
    // expression of the missed-cleavage count and the sequence strings made this a serial step of 20 - 40 ms per file.
    // The distinct peptides and proteins are counted in parallel, too: each range notes its sequence strings (in one
    // buffer) and accessions in buckets by hash, and the buckets are counted one per thread. Merging sets of strings
    // per range and freeing them afterwards took 0.6 s on one thread for the 250,000 top hits of a timsTOF run.
    //
    // The mass errors of the tolerance estimate (with_errors), of the top hits of the target PSMs (a decoy is a wrong
    // match: its errors are noise), signed (observed minus theoretical, ppm): of the precursor, corrected by the isotope
    // error of the hit (as the search matched it; @p isotope_observed_minus_theoretical: the convention of the reported
    // 'isotope_error'), without the hits with a delta mass of an open search; and of the annotated b and y fragment peaks.
    TopHitSummary summarizeTopHits(const PeptideIdentificationList& peptide_ids, const std::string& enzyme, bool with_errors,
                                   bool isotope_observed_minus_theoretical = true)
    {
      const DigestionEnzymeProtein* protease = ProteaseDB::getInstance()->getEnzyme(enzyme);
      const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      const UInt target_decoy_index = registry.getIndex(Constants::UserParam::TARGET_DECOY); // UInt(-1): exists nowhere
      const UInt delta_mass_index = registry.getIndex("DeltaMass");
      const UInt isotope_error_index = registry.getIndex(Constants::UserParam::ISOTOPE_ERROR);
      static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
      const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
      // Declared without OpenMP, too: with the -fopenmp-simd fallback the compiler still parses the clauses below.
      int num_parts = 1;
#ifdef _OPENMP
      num_parts = static_cast<int>(std::clamp<SignedSize>(n / 512, 1, omp_get_max_threads()));
#endif
      const Size num_buckets = static_cast<Size>(num_parts);
      struct Part
      {
        TopHitSummary summary;
        std::string sequences;                                             // the sequence strings, one after the other
        std::vector<std::vector<std::pair<Size, Size>>> peptide_buckets;   // (offset, length) in sequences, by hash
        std::vector<std::vector<const std::string*>> protein_buckets;      // accessions of the peptide evidences, by hash
      };
      std::vector<Part> parts(num_parts);
      std::vector<std::exception_ptr> errors(num_parts); // an exception must not leave the parallel region
#pragma omp parallel for schedule(static, 1) num_threads(num_parts)
      for (int part = 0; part < num_parts; ++part)
      try
      {
        Part& p = parts[part];
        TopHitSummary& s = p.summary;
        p.peptide_buckets.resize(num_buckets);
        p.protein_buckets.resize(num_buckets);
        const std::hash<std::string_view> hash;
        EnzymaticDigestion digestor;
        digestor.setEnzyme(protease);
        std::unordered_map<std::string, Size> missed_cleavages_of; // by unmodified sequence
        std::vector<double> residue_masses; // of a hit, cumulative (prefix sums, for the b and y ions)
        for (SignedSize i = n * part / num_parts; i < n * (part + 1) / num_parts; ++i)
        {
          const PeptideIdentification& pid = peptide_ids[i];
          if (pid.getHits().empty()) { continue; }
          ++s.matched;
          const PeptideHit& top = pid.getHits().front();
          const std::string sequence = top.getSequence().toString();
          p.peptide_buckets[hash(sequence) % num_buckets].emplace_back(p.sequences.size(), sequence.size());
          p.sequences += sequence;
          for (const auto& ev : top.getPeptideEvidences())
          {
            const std::string& accession = ev.getProteinAccession();
            p.protein_buckets[hash(accession) % num_buckets].push_back(&accession);
          }
          if (top.getCharge() > 0) { ++s.charges[top.getCharge()]; }
          // the target/decoy type as PeptideHit::getTargetDecoyType() (by the index of the meta value; case-insensitive;
          // another value is neither, where getTargetDecoyType() throws)
          std::string td;
          if (top.metaValueExists(target_decoy_index)) { td = StringUtils::toLowered(top.getMetaValue(target_decoy_index).toString()); }
          const bool is_target = td == "target" || td == "target+decoy";
          if (with_errors)
          {
            // Precursor error: skip hits with unresolved charge (0) - getMZ() would throw - and, in an open search, those
            // with a delta mass (a modification the sequence lacks, not a mass error; as OpenSearchModificationAnalysis
            // tells unmodified PSMs).
            const bool shifted = top.metaValueExists(delta_mass_index)
                                 && !OpenSearchModificationAnalysis::isUnmodified(static_cast<double>(top.getMetaValue(delta_mass_index)));
            if (is_target && !shifted && top.getCharge() > 0)
            {
              // the isotope offset the search matched at (theoretical minus observed)
              int offset = 0;
              if (top.metaValueExists(isotope_error_index))
              {
                const int reported = static_cast<int>(top.getMetaValue(isotope_error_index));
                offset = isotope_observed_minus_theoretical ? -reported : reported;
              }
              const double corrected_mz = monoisotopicPrecursorMZ(pid.getMZ(), offset, static_cast<double>(top.getCharge()));
              s.precursor_errors.push_back(Math::getPPM(corrected_mz, top.getSequence().getMZ(top.getCharge())));
            }
            // Fragment errors: of the annotated b and y peaks ("b3+", "y12++"; not those with a neutral loss or of other
            // ions), against their m/z from the hit's sequence: (prefix or suffix residue masses + terminal modification
            // [+ H2O] + charge protons) / charge, the masses TheoreticalSpectrumGenerator gives the annotated ions
            if (is_target && !top.getPeakAnnotations().empty())
            {
              const AASequence& sequence = top.getSequence();
              const Size length = sequence.size();
              residue_masses.assign(length + 1, 0.0);
              for (Size k = 0; k < length; ++k)
              {
                residue_masses[k + 1] = residue_masses[k] + sequence[k].getMonoWeight(Residue::Internal);
              }
              const double n_term = sequence.hasNTerminalModification() ? sequence.getNTerminalModification()->getDiffMonoMass() : 0.0;
              const double c_term = sequence.hasCTerminalModification() ? sequence.getCTerminalModification()->getDiffMonoMass() : 0.0;
              bool any = false;
              for (const PeptideHit::PeakAnnotation& pa : top.getPeakAnnotations())
              {
                const std::string& name = pa.annotation;
                bool prefix = false;
                Size ordinal = 0;
                if (name.empty() || (name[0] != 'b' && name[0] != 'y') || pa.charge < 1
                    || !FragmentIonLikelihoodModel::parseIonName(name, prefix, ordinal) || ordinal >= length
                    || name.find_first_not_of('+', name.find_first_not_of("0123456789", 1)) != std::string::npos)
                {
                  continue; // another ion, or one with a neutral loss ("y3-H2O1+")
                }
                const double neutral = prefix ? residue_masses[ordinal] + n_term
                                              : residue_masses[length] - residue_masses[length - ordinal] + c_term + water_mass;
                const double theoretical_mz = (neutral + pa.charge * Constants::PROTON_MASS_U) / pa.charge;
                s.fragment_errors.push_back(Math::getPPM(pa.mz, theoretical_mz));
                any = true;
              }
              if (any) { ++s.fragment_errors_psms; }
            }
          }
          // counted once per sequence of the range: the regular expression takes a block from boost's process-wide
          // memory cache per search, which the threads contended for (93 ms at 16 and 64 threads for a timsTOF run)
          std::string unmodified = top.getSequence().toUnmodifiedString();
          auto cleavages = missed_cleavages_of.find(unmodified);
          if (cleavages == missed_cleavages_of.end())
          {
            const Size count = digestor.countInternalCleavageSites(unmodified);
            cleavages = missed_cleavages_of.emplace(std::move(unmodified), count).first;
          }
          ++s.missed_cleavages[cleavages->second];
          if (td == "decoy") { ++s.n_decoy; }
          else if (is_target) { ++s.n_target; }
        }
      }
      catch (...)
      {
        errors[part] = std::current_exception();
      }
      for (const std::exception_ptr& error : errors)
      {
        if (error) { std::rethrow_exception(error); }
      }

      // The distinct strings of each bucket (equal strings have equal hashes, so they share a bucket), and the buffers of
      // each range freed by a thread of the team
      std::vector<Size> unique_peptides(num_buckets, 0), unique_proteins(num_buckets, 0);
      std::vector<std::exception_ptr> bucket_errors(num_buckets);
#pragma omp parallel num_threads(num_parts)
      {
        std::vector<std::string_view> strings;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(num_buckets); ++sb)
        try
        {
          const Size b = static_cast<Size>(sb);
          strings.clear();
          for (const Part& p : parts)
          {
            for (const auto& [offset, length] : p.peptide_buckets[b]) { strings.emplace_back(p.sequences.data() + offset, length); }
          }
          std::sort(strings.begin(), strings.end());
          unique_peptides[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
          strings.clear();
          for (const Part& p : parts)
          {
            for (const std::string* accession : p.protein_buckets[b]) { strings.emplace_back(*accession); }
          }
          std::sort(strings.begin(), strings.end());
          unique_proteins[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
        }
        catch (...)
        {
          bucket_errors[sb] = std::current_exception();
        }
#pragma omp for schedule(static, 1)
        for (int part = 0; part < num_parts; ++part)
        {
          std::string().swap(parts[part].sequences);
          decltype(Part::peptide_buckets)().swap(parts[part].peptide_buckets);
          decltype(Part::protein_buckets)().swap(parts[part].protein_buckets);
        }
      }

      for (const std::exception_ptr& error : bucket_errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
      Size num_precursor_errors = 0, num_fragment_errors = 0;
      for (const Part& p : parts)
      {
        num_precursor_errors += p.summary.precursor_errors.size();
        num_fragment_errors += p.summary.fragment_errors.size();
      }
      TopHitSummary all = std::move(parts.front().summary);
      all.unique_peptides = std::accumulate(unique_peptides.begin(), unique_peptides.end(), Size(0));
      all.unique_proteins = std::accumulate(unique_proteins.begin(), unique_proteins.end(), Size(0));
      all.precursor_errors.reserve(num_precursor_errors);
      all.fragment_errors.reserve(num_fragment_errors);
      for (int part = 1; part < num_parts; ++part)
      {
        TopHitSummary& s = parts[part].summary;
        all.matched += s.matched;
        all.n_target += s.n_target;
        all.n_decoy += s.n_decoy;
        all.fragment_errors_psms += s.fragment_errors_psms;
        for (const auto& [charge, count] : s.charges) { all.charges[charge] += count; }
        for (const auto& [cleavages, count] : s.missed_cleavages) { all.missed_cleavages[cleavages] += count; }
        all.precursor_errors.insert(all.precursor_errors.end(), s.precursor_errors.begin(), s.precursor_errors.end());
        all.fragment_errors.insert(all.fragment_errors.end(), s.fragment_errors.begin(), s.fragment_errors.end());
      }
      return all;
    }
    // Math::median() and Math::MAD() of @p values (reordered), by selection instead of sorts: the same order statistics,
    // so the same values (the sorts of the 250,000 errors of a timsTOF run took 20-30 ms each). Unless a value is not
    // finite (a NaN is not ordered; an infinite median gives NaN deviations) or a negative zero (equal to 0, but written
    // differently), where the result may depend on the algorithm or on the order the sort left: then Math itself.
    std::pair<double, double> medianAndMADBySelection(std::vector<double>& values)
    {
      const auto ambiguous = [](double v) { return !std::isfinite(v) || (v == 0.0 && std::signbit(v)); };
      if (values.empty() || std::any_of(values.begin(), values.end(), ambiguous))
      {
        const double median = Math::median(values.begin(), values.end());
        return {median, Math::MAD(values.begin(), values.end(), median)};
      }
      // the median of finite values without a negative zero; their deviations are not negative zeros and not NaN
      // Large lists (the b/y fragment errors of a timsTOF run: millions) with the OpenMP threads: the selection by
      // libstdc++'s parallel nth_element (the same k-th smallest value), the deviations in parallel.
      constexpr Size parallel_from = 65536;
      const auto median_of = [](std::vector<double>& v)
      {
        const auto upper = v.begin() + static_cast<SignedSize>(v.size() / 2);
#ifdef PROSE_GNU_PARALLEL
        if (v.size() >= parallel_from && omp_get_max_threads() > 1) { __gnu_parallel::nth_element(v.begin(), upper, v.end()); }
        else
#endif
        { std::nth_element(v.begin(), upper, v.end()); }
        if (v.size() % 2 == 1) { return *upper; }
        return (*std::max_element(v.begin(), upper) + *upper) / 2.0;
      };
      const double median = median_of(values);
      std::vector<double> diffs(values.size());
#pragma omp parallel for schedule(static) if (values.size() >= parallel_from)
      for (SignedSize i = 0; i < static_cast<SignedSize>(values.size()); ++i) { diffs[i] = fabs(values[i] - median); }
      return {median, median_of(diffs)};
    }

  // =====================================================================
  // Recompute result-level stats from a FINAL (post-rescoring/post-FDR) PSM list.
  // =====================================================================
    // Spectra with at least one PSM. The PSMs of co-isolated precursors (chimeric:enabled) belong to the spectrum
    // (scan_index) of the selected precursor, which counts once.
    Size countMatchedSpectra(const PeptideIdentificationList& peptide_ids)
    {
      Size matched = 0;
      std::set<unsigned int> selected_scans, coisolated_scans;
      for (const PeptideIdentification& pid : peptide_ids)
      {
        if (pid.getHits().empty()) continue;
        const PeptideHit& top = pid.getHits().front();
        const bool coisolated = top.metaValueExists(ProSEAlgorithm::CHIMERIC_RANK)
                                && static_cast<int>(top.getMetaValue(ProSEAlgorithm::CHIMERIC_RANK)) == 2;
        if (!pid.metaValueExists("scan_index"))
        {
          if (!coisolated) ++matched;
          continue;
        }
        const unsigned int scan = static_cast<unsigned int>(pid.getMetaValue("scan_index"));
        (coisolated ? coisolated_scans : selected_scans).insert(scan);
        if (!coisolated) ++matched;
      }
      for (unsigned int scan : coisolated_scans)
      {
        if (selected_scans.count(scan) == 0) ++matched;
      }
      return matched;
    }

    // countMatchedSpectra(), given @p with_hits, the number of identifications with hits: the same unless a top hit is
    // of a co-isolated precursor (chimeric:enabled), which is checked first (in parallel: no set of scans per search)
    Size matchedSpectra(const PeptideIdentificationList& peptide_ids, Size with_hits)
    {
      const UInt rank_index = MetaInfoInterface::metaRegistry().getIndex(ProSEAlgorithm::CHIMERIC_RANK);
      if (rank_index == UInt(-1)) { return with_hits; } // (a name that was never registered is on no hit)
      bool any_coisolated = false;
#pragma omp parallel for schedule(static) reduction(||: any_coisolated)
      for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
      {
        const auto& hits = peptide_ids[i].getHits();
        if (!hits.empty() && hits.front().metaValueExists(rank_index)
            && static_cast<int>(hits.front().getMetaValue(rank_index)) == 2) { any_coisolated = true; }
      }
      return any_coisolated ? countMatchedSpectra(peptide_ids) : with_hits;
    }
  }

  void ProSEAlgorithm::updateFinalStats(RunStatistics& stats,
                                        const PeptideIdentificationList& peptide_ids,
                                        const std::string& enzyme,
                                        bool fdr_applied)
  {
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    const TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme, false);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    stats.charge_histogram = summary.charges;
    stats.missed_cleavage_histogram = summary.missed_cleavages;
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    stats.fdr_applied = fdr_applied;
    stats.achieved_psm_fdr = fdr_applied ? maxRetainedScore_(peptide_ids) : -1.0;
  }

  void ProSEAlgorithm::releaseInParallel(SearchResult& result)
  {
    for (auto& protein_id : result.protein_ids) { freeInParallel(protein_id.getHits()); }
    std::vector<ProteinIdentification>().swap(result.protein_ids);
    freeInParallel(result.peptide_ids);
    freeInParallel(result.pre_fdr_peptide_ids);
  }

  // =====================================================================
  // Helper: maximum top-hit score among retained PSMs (== achieved FDR after
  // q-value filtering). Returns -1.0 for an empty set.
  // =====================================================================
  double ProSEAlgorithm::maxRetainedScore_(const PeptideIdentificationList& peptide_ids)
  {
    double max_score = -1.0;
    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      max_score = std::max(max_score, pid.getHits().front().getScore());
    }
    return max_score;
  }

  // =====================================================================
  // Helper: fill per-run identification statistics (silent — no logging).
  // Computes target/decoy counts from the (post-FDR) hits it is given, so
  // SearchResult::stats is consistent for any caller; leaves achieved FDR and
  // timing fields untouched (captured at well-defined points in the search paths).
  // =====================================================================
  ProSEAlgorithm::SpectraCounts_ ProSEAlgorithm::countSpectra_(const PeakMap& spectra) const
  {
    // in parallel: one pass over the spectra of a timsTOF run (600,000) took 0.13 s on one thread (memory bound)
    const bool search_outside = param_.getValue("precursor:recorded_charge_outside_range").toString() == "search";
    const int min_charge = static_cast<int>(precursor_min_charge_), max_charge = static_cast<int>(precursor_max_charge_);
    Size ms2 = 0, profile = 0, out_of_range = 0;
#pragma omp parallel reduction(+: ms2, profile, out_of_range) if (spectra.size() >= 16384)
    {
      std::vector<uint16_t> charges;
#pragma omp for schedule(static)
      for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
      {
        const MSSpectrum& s = spectra[i];
        if (s.getMSLevel() != 2) continue;
        ++ms2;
        if (s.getType() == SpectrumSettings::SpectrumType::PROFILE) ++profile;
        if (s.getPrecursors().empty()) continue;
        FragmentIndex::precursorCharges(s.getPrecursors()[0], min_charge, max_charge, charges, search_outside);
        if (charges.empty()) ++out_of_range;
      }
    }
    SpectraCounts_ counts;
    counts.ms2 = ms2;
    counts.profile = profile;
    counts.charge_out_of_range = out_of_range;
    return counts;
  }

  void ProSEAlgorithm::preprocessSearchedSpectra_(PeakMap& spectra, const std::string& input,
                                                   FragmentIonLikelihoodModel::PeakLists& ion_evidence,
                                                   PeakMap* evidence_spectra, PeakMap* query_spectra,
                                                   PeakCharges_* charges) const
  {
    // the charges only for a 'curve' recalibration, which corrects the peaks with them (applyRecalibration_())
    if (charges != nullptr) { *charges = PeakCharges_(); }
    if (!(calibration_enabled_ && calibration_recalibrate_ == "curve" && query_spectra != nullptr && fragment_mass_tolerance_unit_ == "ppm"))
    {
      charges = nullptr;
    }
    const SpectraCounts_ counts = countSpectra_(spectra);
    const std::string of = input.empty() ? std::string() : " of '" + input + "'";
    if (counts.ms2 == 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << (input.empty() ? std::string("the input") : "'" + input + "'")
                      << " holds no MS2 spectra (e.g. an MS1-only file): "
                      << "nothing is searched and no peptides are identified." << std::endl;
    }
    if (counts.profile > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.profile << " of " << counts.ms2 << " MS2 spectra" << of
                      << " are recorded as profile spectra. ProSE expects centroided spectra and searches profile "
                      << "spectra as they are, which identifies few peptides, if any. Centroid them first (e.g. with "
                      << "PeakPickerHiRes, or peak picking when converting the raw files)." << std::endl;
    }
    if (counts.charge_out_of_range > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.charge_out_of_range << " of " << counts.ms2
                      << " MS2 spectra" << of << " have a precursor charge outside of precursor:min_charge.."
                      << "max_charge (" << precursor_min_charge_ << ".." << precursor_max_charge_
                      << ") and are not searched." << std::endl;
    }
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm", deisotope_requested_,
                       peaks_keep_n_, peaks_window_top_, peaks_window_type_, deisotoping_,
                       self_trained_ion_priors_ ? &ion_evidence : nullptr, ion_prior_scored_peaks_, evidence_spectra,
                       query_spectra, charges);
  }

  void ProSEAlgorithm::collectRunStatistics_(
      const PeakMap& spectra,
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids,
      RunStatistics& stats) const
  {
    const SpectraCounts_ counts = countSpectra_(spectra);
    stats.ms2_spectra = counts.ms2;
    stats.profile_spectra = counts.profile;
    stats.charge_out_of_range_spectra = counts.charge_out_of_range;
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    // Target/decoy of the FINAL (post-FDR) hits, so SearchResult::stats stays consistent
    // for library callers that read it without calling updateFinalStats().
    // Mass errors for the tolerance estimate: see summarizeTopHits() (signed, isotope-corrected, of the target PSMs).
    TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme_, true, isotope_error_observed_minus_theoretical_);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    for (const auto& [charge, count] : summary.charges) { stats.charge_histogram[charge] += count; }
    for (const auto& [cleavages, count] : summary.missed_cleavages) { stats.missed_cleavage_histogram[cleavages] += count; }
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    vector<double>& precursor_errors = summary.precursor_errors;
    vector<double>& fragment_errors = summary.fragment_errors;
    const Size fragment_errors_psms = summary.fragment_errors_psms;

    // -- Per-run tolerance estimate: the median of the signed errors (the bias) and their MAD (the spread; unscaled).
    //    The recommended (symmetric) tolerance covers the bias plus three standard deviations of a normal distribution,
    //    estimated robustly as 1.4826 * MAD: ceil(|median| + 3 * 1.4826 * MAD). Of the fragments it is an estimate
    //    within the searched tolerance: peaks beyond it were not matched. --
    const Size min_psms_for_estimation = 10;
    const auto recommended_tolerance = [](double median, double mad) { return std::ceil(std::fabs(median) + 3.0 * 1.4826 * mad); };
    if (precursor_errors.size() >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(precursor_errors); // = Math::median() / Math::MAD()
      stats.prec_err_median = med;
      stats.prec_err_mad = mad;
      stats.prec_err_recommended = recommended_tolerance(med, mad);
      stats.prec_tol_valid = true;
    }
    // With a 'curve' model the fragment annotations report the m/z as searched (corrected; see
    // appliedFragmentRecalibration()), so the fragment errors are those after the curve
    stats.frag_err_after_curve = false;
    if (!protein_ids.empty())
    {
      const ProteinIdentification::SearchParameters& sp = protein_ids[0].getSearchParameters();
      stats.frag_err_after_curve = sp.getMetaValue("calibration:recalibrate", DataValue("")).toString() == "curve"
                                   && sp.getMetaValue("calibration:fragment_model", DataValue("")).toString().find("cubic") != std::string::npos;
    }
    if (fragment_errors_psms >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(fragment_errors);
      stats.frag_err_median = med;
      stats.frag_err_mad = mad;
      stats.frag_err_recommended = recommended_tolerance(med, mad);
      stats.frag_tol_valid = true;
    }

    // -- calibration: the records of recordRecalibration_() (measured errors; recalibration), for the run summary --
    stats.recalibration.clear();
    if (!protein_ids.empty() && protein_ids[0].getSearchParameters().metaValueExists("calibration:precursor_error_offset"))
    {
      const ProteinIdentification::SearchParameters& sp = protein_ids[0].getSearchParameters();
      std::vector<std::string> keys;
      sp.getKeys(keys);
      for (const std::string& key : keys)
      {
        if (key.compare(0, 12, "calibration:") == 0) { stats.recalibration.emplace_back(key.substr(12), sp.getMetaValue(key)); }
      }
    }
  }

  // =====================================================================
  // Self-trained ion priors: learn presence, intensity rank and mass error
  // of fragment ions from the confident PSMs of this run, one model per
  // half of the spectra, and annotate every hit with the other half's model.
  // =====================================================================
  // static
  AASequence ProSEAlgorithm::reversedNoiseSequence_(const AASequence& sequence)
  {
    // Reverse all but the C-terminal residue and keep every residue's modification, so the noise
    // hypothesis shares composition, mass and enzymatic C-terminus with the peptide.
    const Size n = sequence.size();
    if (n < 3) return sequence;
    // The residues are those of the sequence: fromResidues() takes them as they are, where operator+=(const Residue*)
    // looks every residue up in the ResidueDB under a critical section (the threads of annotateIonPriors_() spent
    // a third of their time in that section in stack samples of a 27-run search).
    std::vector<const Residue*> residues;
    residues.reserve(n);
    for (Size i = n - 1; i-- > 0;) { residues.push_back(&sequence[i]); }
    residues.push_back(&sequence[n - 1]);
    AASequence reversed = AASequence::fromResidues(std::move(residues));
    if (sequence.hasNTerminalModification()) { reversed.setNTerminalModification(sequence.getNTerminalModification()); }
    if (sequence.hasCTerminalModification()) { reversed.setCTerminalModification(sequence.getCTerminalModification()); }
    return reversed;
  }

  namespace
  {
  // =====================================================================
  // Helper: run calibration pass on a subset of spectra
  // =====================================================================
    /// The signed errors (ppm, or Da) of the peaks of @p spec closest to the theoretical peaks of @p theo within
    /// @p tolerance (unmatched theoretical peaks contribute nothing); both spectra sorted by m/z. With @p mzs, the
    /// theoretical m/z of each matched peak, aligned with @p errors. With @p all, every peak within @p tolerance of a
    /// theoretical peak instead of the closest one (calibration:recalibrate curve).
    void matchedFragmentErrors(const MSSpectrum& spec, const PeakSpectrum& theo, double tolerance, bool ppm, std::vector<float>& errors,
                               std::vector<float>* mzs = nullptr, bool all = false)
    {
      for (const Peak1D& t : theo)
      {
        const double tol = ppm ? t.getMZ() * tolerance * 1e-6 : tolerance;
        if (all)
        {
          for (auto it = spec.MZBegin(t.getMZ() - tol); it != spec.end() && it->getMZ() <= t.getMZ() + tol; ++it)
          {
            const double e = it->getMZ() - t.getMZ();
            errors.push_back(static_cast<float>(ppm ? e / t.getMZ() * 1e6 : e));
            if (mzs != nullptr) mzs->push_back(static_cast<float>(t.getMZ()));
          }
          continue;
        }
        const Int i = spec.findNearest(t.getMZ(), tol);
        if (i < 0) continue;
        const double e = spec[i].getMZ() - t.getMZ();
        errors.push_back(static_cast<float>(ppm ? e / t.getMZ() * 1e6 : e));
        if (mzs != nullptr) mzs->push_back(static_cast<float>(t.getMZ()));
      }
    }

    /// A cubic polynomial in t, the position of s in [lo, hi] scaled to [-1, 1] (the abscissae of its fit)
    struct Cubic
    {
      double c[4]{0.0, 0.0, 0.0, 0.0};
      double lo{0.0}, hi{0.0};
      double operator()(double s) const
      {
        const double t = 2.0 * (s - lo) / (hi - lo) - 1.0;
        return c[0] + t * (c[1] + t * (c[2] + t * c[3]));
      }
    };

    /// The least-squares cubic through the points (s, y); false when it cannot be fitted
    bool fitCubic(const std::vector<double>& s, const std::vector<double>& y, Cubic& out)
    {
      if (s.size() < 4) return false;
      const auto [lo, hi] = std::minmax_element(s.begin(), s.end());
      if (!(*hi > *lo)) return false;
      Cubic f;
      f.lo = *lo;
      f.hi = *hi;
      Eigen::Matrix4d a = Eigen::Matrix4d::Zero();
      Eigen::Vector4d b = Eigen::Vector4d::Zero();
      for (Size i = 0; i < s.size(); ++i)
      {
        const double t = 2.0 * (s[i] - f.lo) / (f.hi - f.lo) - 1.0;
        const Eigen::Vector4d p(1.0, t, t * t, t * t * t);
        a += p * p.transpose();
        b += p * y[i];
      }
      const Eigen::LDLT<Eigen::Matrix4d> ldlt(a);
      if (ldlt.info() != Eigen::Success || !ldlt.isPositive() || ldlt.rcond() < 1e-12) return false; // (m/z at fewer than 4 places)
      const Eigen::Vector4d c = ldlt.solve(b);
      if (!c.allFinite()) return false;
      for (int k = 0; k < 4; ++k) f.c[k] = c[k];
      out = f;
      return true;
    }

    /// The error curve of the fragment peaks of one RT chunk (calibration:recalibrate curve). @p errors (ppm) holds
    /// every peak within the collection window of a theoretical fragment at @p mzs, most of them random matches: the
    /// mode of the errors in up to 20 m/z bins gives a first cubic in sqrt(m/z), which is refitted three times to the
    /// peaks near it. False when there are too few peaks; @p points: the peaks of the last fit.
    bool fitErrorCurve(const std::vector<double>& errors, const std::vector<double>& mzs, Cubic& curve, Size& points)
    {
      const Size n = errors.size();
      const Size n_bins = std::min<Size>(20, n / 30);
      if (n_bins < 4) return false;
      std::vector<Size> order(n);
      std::iota(order.begin(), order.end(), Size(0));
      std::sort(order.begin(), order.end(), [&mzs](Size x, Size y) { return mzs[x] < mzs[y]; });
      if (!(mzs[order.front()] > 0.0)) return false; // (sqrt)
      std::vector<double> bin_s, bin_mode;
      for (Size bin = 0; bin < n_bins; ++bin)
      {
        const Size begin = n * bin / n_bins, end = n * (bin + 1) / n_bins;
        double lo = errors[order[begin]], hi = lo;
        for (Size i = begin; i < end; ++i)
        {
          lo = std::min(lo, errors[order[i]]);
          hi = std::max(hi, errors[order[i]]);
        }
        lo = std::floor(lo);
        if (!(hi - lo < 1e5)) return false; // (not the errors of a ppm window)
        std::vector<Size> hist(static_cast<Size>(hi - lo) + 1, 0); // 1 ppm per bin
        for (Size i = begin; i < end; ++i) ++hist[static_cast<Size>(errors[order[i]] - lo)];
        Size best = 0, best_count = 0;
        for (Size k = 0; k < hist.size(); ++k) // the (first) maximum of the histogram summed over 5 ppm
        {
          Size count = 0;
          for (Size j = (k < 2 ? 0 : k - 2); j < std::min(hist.size(), k + 3); ++j) count += hist[j];
          if (count > best_count) { best_count = count; best = k; }
        }
        bin_s.push_back(std::sqrt(mzs[order[(begin + end) / 2]]));
        bin_mode.push_back(lo + static_cast<double>(best) + 0.5);
      }
      Cubic f;
      if (!fitCubic(bin_s, bin_mode, f)) return false;
      {
        // a bin whose mode is a cluster of random matches lies far from the curve of the others: fitted once more without
        // ponytail: one round, and a third of the bins off the curve still bends it; RANSAC over the bins if that is seen
        std::vector<double> off(n_bins);
        for (Size k = 0; k < n_bins; ++k) off[k] = std::fabs(bin_mode[k] - f(bin_s[k]));
        std::vector<double> sorted_off = off;
        const double limit = std::max(15.0, 3.0 * 1.4826 * medianAndMADBySelection(sorted_off).first);
        std::vector<double> kept_s, kept_mode;
        for (Size k = 0; k < n_bins; ++k)
        {
          if (off[k] > limit) continue;
          kept_s.push_back(bin_s[k]);
          kept_mode.push_back(bin_mode[k]);
        }
        if (kept_s.size() < n_bins && kept_s.size() >= 4) fitCubic(kept_s, kept_mode, f); // (f stays if this fails)
      }
      // ponytail: the cut follows the spread of the peaks near the curve (3 robust sigma, 2 to 15 ppm), starting at
      // 15 ppm; no model of the random matches, which lie on both sides of the curve
      double cut = 15.0;
      std::vector<double> s, e, residuals;
      for (int round = 0; round < 3; ++round)
      {
        s.clear();
        e.clear();
        residuals.clear();
        for (Size i = 0; i < n; ++i)
        {
          const double si = std::sqrt(mzs[i]);
          const double r = errors[i] - f(si);
          if (std::fabs(r) > cut) continue;
          s.push_back(si);
          e.push_back(errors[i]);
          residuals.push_back(r);
        }
        if (!fitCubic(s, e, f)) return false;
        for (Size i = 0; i < s.size(); ++i) residuals[i] = e[i] - f(s[i]); // about the new curve
        cut = std::clamp(3.0 * 1.4826 * medianAndMADBySelection(residuals).second, 2.0, 15.0);
      }
      curve = f;
      points = s.size();
      return true;
    }
  }

  ProSEAlgorithm::CalibrationResult_
  ProSEAlgorithm::runCalibrationPass_(
      PeakMap& spectra,
      FragmentIndex& fragment_index,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const PeakMap* query_spectra,
      const std::string& decoy_string,
      bool decoy_is_prefix,
      const PeakCharges_* charges) const
  {
    CalibrationResult_ result;

    const bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    const bool precursor_ppm = (precursor_mass_tolerance_unit_ == "ppm");

    // The pass searches with windows wide enough to measure the error (calibrationPassWindow_(); the windows of the main
    // search are derived from the measurement below and never exceed the configured ones). Query-time parameters of the
    // index only, restored on exit.
    const double pass_lower = calibrationPassWindow_(precursor_mass_tolerance_lower_, precursor_ppm);
    const double pass_upper = calibrationPassWindow_(precursor_mass_tolerance_upper_, precursor_ppm);
    const double pass_fragment = calibrationPassWindow_(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
    // calibration:recalibrate curve: an error that varies with m/z by more than the pass window is cut off by it, so
    // every peak within 150 ppm of a fragment is collected and the peaks near the fitted curve are kept (below)
    const bool curve = calibration_recalibrate_ == "curve" && fragment_mass_tolerance_unit_ppm && query_spectra != nullptr;
    if (calibration_recalibrate_ == "curve" && !curve)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: calibration:recalibrate curve needs a ppm fragment tolerance and the peaks as measured "
                         "(fragment:query_spectrum raw); a line in m/z is fitted instead." << std::endl;
    }
    const double collect_fragment = curve ? std::max(pass_fragment, 150.0) : pass_fragment;
    const Param fi_params_saved = fragment_index.getParameters();
    {
      Param p = fi_params_saved;
      p.setValue("precursor:mass_tolerance_lower", pass_lower);
      p.setValue("precursor:mass_tolerance_upper", pass_upper);
      p.setValue("fragment:mass_tolerance", pass_fragment);
      fragment_index.setParameters(p);
    }
    const OnScopeExit restore_fi([&fragment_index, &fi_params_saved]()
    {
      try { fragment_index.setParameters(fi_params_saved); }
      catch (...) {} // only after another exception: the caller restores the parameters of its index
    });

    // The subset: the spectra of the highest TIC (the best-measured ones), taken per tenth of the run (the spectra are in
    // acquisition order) so that a mass offset that changes over the run is represented in the sample
    vector<pair<double, Size>> tic_index(spectra.size());
#pragma omp parallel for schedule(static) if (spectra.size() >= 4096)
    for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
    {
      double tic = 0;
      for (const auto& p : spectra[i]) { tic += p.getIntensity(); }
      tic_index[i] = {tic, static_cast<Size>(i)};
    }
    vector<Size> subset;
    {
      const Size blocks = std::min<Size>(10, spectra.size());
      for (Size b = 0; b < blocks; ++b)
      {
        const Size begin = spectra.size() * b / blocks, end = spectra.size() * (b + 1) / blocks;
        const Size take = std::min<Size>(end - begin, std::max<Size>(1, static_cast<Size>(std::llround((end - begin) * calibration_subset_ratio_))));
        std::partial_sort(tic_index.begin() + begin, tic_index.begin() + begin + take, tic_index.begin() + end,
                          [](const pair<double, Size>& a, const pair<double, Size>& c) { return a.first > c.first; });
        for (Size i = begin; i < begin + take; ++i) { subset.push_back(tic_index[i].second); }
      }
    }
    // ('curve': in acquisition order, as the passes of its first version scored it; otherwise as upstream, which
    // only decides between hits of the same score at the cut below)
    if (curve) { std::sort(subset.begin(), subset.end()); }
    const Size subset_size = subset.size();
    // Scored in the order of the precursor mass, as the main search (scoringOrder_()): neighbouring iterations probe the
    // same precursor-window slices of the fragment index, which stay in the caches. The hits are collected in the order
    // of the subset (per tenth of the run, by TIC), as one thread scored them before: the sort by score below is not
    // stable, so this order decides which of two hits with the same score at the cut survives; it depends neither on
    // the threads nor on the chunks.
    std::vector<Size> order; // positions in subset, by precursor mass
    {
      std::vector<SignedSize> position_in_subset(spectra.size(), -1);
      for (Size k = 0; k < subset_size; ++k) position_in_subset[subset[k]] = static_cast<SignedSize>(k);
      order.reserve(subset_size);
      for (SignedSize i : scoringOrder_(spectra, query_spectra))
      {
        if (position_in_subset[i] >= 0) order.push_back(static_cast<Size>(position_in_subset[i]));
      }
    }

    OPENMS_LOG_INFO << "[ProSE] Calibration: scoring " << subset_size << " / " << spectra.size()
                    << " spectra (top TIC per tenth of the run)..." << std::endl;

    // Score subset and collect errors from the best hit per spectrum, with the ion series
    // the fragment index was built from (e.g. c/z+1 for ETD, where b/y would match nothing)
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // Collect per-spectrum best hits with scores and errors
    // (rt, prec_mz [theoretical] and frag_mzs [theoretical, aligned with frag_errors] feed the models of calibration:recalibrate;
    // decoy: the best hit is a peptide of a decoy protein, which marks the RT chunks of the precursor model whose PSMs are
    // mostly false matches)
    struct CalHit { double score; double prec_error; std::vector<float> frag_errors; double rt; double prec_mz; std::vector<float> frag_mzs; bool decoy; };
    // calibration:recalibrate curve scores the subset again as corrected by the last curve until the curve is stable (below)
    Size cal_hits_total = 0;
    vector<double> precursor_errors, precursor_mzs, precursor_rts;
    vector<char> precursor_decoys;
    vector<double> fragment_errors, fragment_mzs, fragment_rts;
    MzRecalibration_ fragment_curve, previous_curve;
    // 'curve': the first pass only has to give the curve that the next pass corrects the subset with, so it scores an
    // evenly spaced sample of about 10000 spectra of the subset; again all of them if that sample fails. Each later pass
    // scores the whole subset as corrected by the last curve and fits the curve again (to the errors as measured), until
    // it moves by at most curve_tolerance ppm where the fragments are, or after max_passes.
    // ponytail: fixed sample size, tolerance and pass limit; parameters if runs need others
    constexpr Size first_pass_sample = 10000;
    constexpr double curve_tolerance = 2.0; // ppm (the residual spread of a timsTOF run is about 6 ppm)
    constexpr int max_passes = 3;
    Size sample = curve && subset_size > first_pass_sample ? first_pass_sample : 0; // (0: all)
    // The later passes score each spectrum of the subset as corrected by the last curve (correcting; inactive in the first
    // pass): the searched peaks at the m/z they were measured at (with the charges of the deisotoper) and the peaks as
    // measured, in copies of this thread; the spectra stay as measured, and the errors are read from them.
    const bool later_passes = curve && charges != nullptr && charges->searched.size() == spectra.size();
    MzRecalibration_ correcting;
    for (int pass = 0; ; ++pass)
    {
      const Size pass_sample = pass == 0 ? sample : 0;
      // At most one hit per spectrum of the subset, in its slot (written by one iteration only): collected in subset order
      // after the loop
      vector<CalHit> slots(subset_size);
      std::vector<unsigned char> filled(subset_size, 0);

      // Parallelize over the calibration subset, mirroring the main scoring loop
      // (scoreSpectraAgainstIndex_). Each iteration is independent: querySpectrum and
      // the shared TheoreticalSpectrumGenerators expose const, thread-safe methods (the
      // main loop already calls them concurrently), and every working variable below is
      // loop-local. The only cross-thread write is the push into cal_hits, guarded by a
      // critical section. Without this the calibration pass ran single-threaded: on a
      // many-core machine the 10% subset took about as long as the entire parallel
      // search, roughly doubling wall time for no extra useful work. Downstream is
      // order-independent — cal_hits is sorted by score and the error vectors are sorted
      // before quantiles — so parallel insertion order does not change the result.
      // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
      // leave the parallel region, which would terminate the process: the first one is rethrown after the loop.
      std::exception_ptr calibration_error;
      std::atomic<bool> calibration_failed{false};
      // Runs of neighbouring masses per core, at least 8 chunks per thread; one spectrum per chunk in peptide-mass mode
      // (as the main loop)
#ifdef _OPENMP
      const int cal_chunk = fragment_index.isPeptideMassMode() ? 1 : static_cast<int>(std::clamp<SignedSize>(
        static_cast<SignedSize>(subset_size) / (SignedSize(8) * std::max(1, omp_get_max_threads())), 1, 64));
#else
      [[maybe_unused]] const int cal_chunk = 1; // (only the OpenMP directive uses it)
#endif
#pragma omp parallel for schedule(dynamic, cal_chunk)
      for (SignedSize oi = 0; oi < (SignedSize)subset_size; ++oi)
      {
        if (calibration_failed.load(std::memory_order_relaxed)) continue;
        if (pass_sample > 0 && order[oi] * pass_sample / subset_size == (order[oi] + 1) * pass_sample / subset_size) continue;
        try
        {
        const Size si = order[oi];
        const Size scan_idx = subset[si];
        // (the peaks as measured: the errors of every pass are read from them)
        const MSSpectrum& measured = query_spectra != nullptr ? (*query_spectra)[scan_idx] : spectra[scan_idx];
        thread_local MSSpectrum corrected_spectrum, corrected_query; // (a later pass: this thread's copies, reused)
        if (correcting.active())
        {
          const MzRecalibration_::AtRT model = correcting.at(spectra[scan_idx].getRT());
          corrected_spectrum = spectra[scan_idx];
          correctAtMeasuredMZ(corrected_spectrum, model, charges->searched[scan_idx]);
          corrected_query = measured;
          correctAtMeasuredMZ(corrected_query, model, {});
        }
        const MSSpectrum& spec = correcting.active() ? corrected_spectrum : spectra[scan_idx];
        const MSSpectrum& query = correcting.active() ? corrected_query : measured;
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);

        // The per-thread containers of the main loop (scoreSpectraAgainstIndex_()): the candidates and the theoretical
        // spectrum keep the capacity the thread's largest spectrum needed instead of taking fresh buffers per spectrum
        ScoringScratch_& scratch = scoringScratch_();
        FragmentIndex::SpectrumMatchesTopN& top_sms = scratch.top_sms;
        top_sms.clear();
        fragment_index.querySpectrum(query, db, top_sms, generators.electronIons(spec));

        // Find the best-scoring hit for this spectrum
        double best_score = 0;
        AASequence best_seq;
        int best_isotope_error = 0;
        uint16_t best_charge = 0;
        UInt32 best_protein = 0; // the protein the fragment index records for the peptide of the best hit

        PeakSpectrum& theo = scratch.theo_spectrum;
        std::vector<uint8_t>& theo_ion_series = scratch.theo_ion_series;

        for (const auto& sms : top_sms.hits_)
        {
          const FragmentIndex::Peptide& peptide = fragment_index.getPeptides()[sms.peptide_idx_];
          AASequence seq = fragment_index.isPeptideMassMode()
            ? fragment_index.reconstructRealizedSubSequence(peptide, db, sms.realized_length_, sms.subset_bitmask_)
            : fragment_index.reconstructModifiedSequence(peptide, db);

          // The calibration PSMs are selected with the search's own score (upstream 825c33bb). The mass-accuracy score
          // favours PSMs with small fragment errors, so it narrows the fragment tolerance estimated from them below
          // (2-10% narrower than with HyperScore on 12 of 14 ppm runs); selecting by HyperScore instead did not
          // change the yield measurably.
          // The peaks of getSpectrum() with the ion series per peak instead of the ion names (no strings, no data
          // arrays), as the main loop: the same scores. The string path only when the generator cannot (losses, ...).
          HyperScore::PSMDetail detail;
          double score;
          if (tsg.getIonSeriesSpectrum(theo, theo_ion_series, seq, 1, scoringMaxCharge_(sms.precursor_charge_)))
          {
            score = mass_accuracy_score_
              ? HyperScore::computeMassAccuracy(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, theo_ion_series, mass_error_sd_ppm_, detail)
              : HyperScore::computeWithDetail(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, theo_ion_series, detail);
          }
          else
          {
            theo.clear(false); // peaks and data arrays; getSpectrum appends, its output is sorted with add_metainfo=true
            theo.getPrecursors().clear(); // (clear(false) keeps the meta data; getSpectrum() appends a precursor per call)
            tsg.getSpectrum(theo, seq, 1, scoringMaxCharge_(sms.precursor_charge_));
            score = mass_accuracy_score_
              ? HyperScore::computeMassAccuracy(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, mass_error_sd_ppm_, detail)
              : HyperScore::computeWithDetail(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, detail);
          }

          if (score > best_score)
          {
            best_score = score;
            best_seq = std::move(seq);
            best_isotope_error = sms.isotope_error_;
            best_charge = sms.precursor_charge_;
            best_protein = peptide.protein_idx;
          }
        }

        if (best_score == 0 || best_seq.empty()) continue;

        // Skip PSMs matched at a non-zero isotope offset. Their precursor m/z carries
        // an extra source of uncertainty (the +1/+2 isotope peak can be ambiguous in
        // the MS1 precursor picking), which inflates the variance of the calibration
        // quantile estimate. iso_err=0 PSMs are the gold-standard "true monoisotopic
        // peak picked" subset — the right population for estimating instrument bias.
        if (best_isotope_error != 0) continue;

        // Compute precursor error (signed), isotope-corrected. FragmentIndex searches
        // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs +
        // isotope_error * C13C12; the observed-to-monoiso m/z correction is
        //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
        // Matches the sign used by postProcessHits_'s PRECURSOR_ERROR_PPM annotation.
        double exp_mz = spec.getPrecursors()[0].getMZ();
        double theo_mz = best_seq.getMZ(best_charge);
        double corrected_exp_mz = monoisotopicPrecursorMZ(exp_mz, best_isotope_error, static_cast<double>(best_charge));
        double prec_err = (precursor_mass_tolerance_unit_ == "ppm")
                            ? Math::getPPM(corrected_exp_mz, theo_mz)
                            : (corrected_exp_mz - theo_mz);

        // The signed errors of the fragment peaks the best hit explains, within the pass window ('curve': of all peaks
        // within the collection window)
        std::vector<float> frag_errors, frag_mzs;
        if (!tsg.getIonSeriesSpectrum(theo, theo_ion_series, best_seq, 1, scoringMaxCharge_(best_charge))) // (the same peaks)
        {
          theo.clear(false);
          theo.getPrecursors().clear();
          tsg.getSpectrum(theo, best_seq, 1, scoringMaxCharge_(best_charge));
        }
        // 'curve' models the error at the m/z a peak was measured at: read it from the peaks as measured where they are
        // at hand (a deisotoped peak has the m/z of the singly charged ion)
        matchedFragmentErrors(curve ? measured : spec, theo, collect_fragment, fragment_mass_tolerance_unit_ppm, frag_errors, &frag_mzs, curve);

        const bool decoy = !decoy_string.empty() && best_protein < db.size()
          && (decoy_is_prefix ? StringUtils::hasPrefix(db[best_protein].identifier, decoy_string)
                              : StringUtils::hasSuffix(db[best_protein].identifier, decoy_string));
        slots[si] = CalHit{best_score, prec_err, std::move(frag_errors), spec.getRT(), theo_mz, std::move(frag_mzs), decoy};
        filled[si] = 1;
        }
        catch (...)
        {
#pragma omp critical (ProSEAlgorithm_calibration_error)
          {
            if (!calibration_error) calibration_error = std::current_exception();
          }
          calibration_failed.store(true, std::memory_order_relaxed);
        }
      }
      if (calibration_error) std::rethrow_exception(calibration_error);
      vector<CalHit> cal_hits;
      cal_hits.reserve(static_cast<Size>(std::count(filled.begin(), filled.end(), 1)));
      for (Size k = 0; k < subset_size; ++k) { if (filled[k]) cal_hits.push_back(std::move(slots[k])); }
      vector<CalHit>().swap(slots);

      // The confident half: the top 50 % by score (random matches inflate the tails), if that leaves calibration_min_psms_
      cal_hits_total = cal_hits.size();
      std::sort(cal_hits.begin(), cal_hits.end(),
                [](const CalHit& a, const CalHit& b) { return a.score > b.score; });
      Size keep = cal_hits.size() / 2;
      if (keep >= calibration_min_psms_ && keep < cal_hits.size())
      {
        cal_hits.resize(keep);
      }

      // The hits the pass window delivers: the signed error e = observed - theoretical lies in [-upper, +lower] (the
      // mass-space window is [observed - lower, observed + upper]; the bounds flip for asymmetric windows)
      for (auto* v : {&precursor_errors, &precursor_mzs, &precursor_rts, &fragment_errors, &fragment_mzs, &fragment_rts}) { v->clear(); }
      precursor_decoys.clear();
      for (const auto& h : cal_hits)
      {
        if (h.prec_error < -pass_upper || h.prec_error > pass_lower) continue; // wrong match
        precursor_errors.push_back(h.prec_error);
        precursor_mzs.push_back(h.prec_mz);
        precursor_rts.push_back(h.rt);
        precursor_decoys.push_back(h.decoy);
        fragment_errors.insert(fragment_errors.end(), h.frag_errors.begin(), h.frag_errors.end());
        fragment_mzs.insert(fragment_mzs.end(), h.frag_mzs.begin(), h.frag_mzs.end());
        fragment_rts.insert(fragment_rts.end(), h.frag_errors.size(), h.rt);
      }

      if (precursor_errors.size() < calibration_min_psms_ && pass_sample > 0)
      {
        sample = 0; // (the sample was too small: all spectra of the subset)
        pass = -1;
        continue;
      }
      if (precursor_errors.size() < calibration_min_psms_)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: insufficient PSMs (" << precursor_errors.size()
                        << " < " << calibration_min_psms_
                        << "), using configured tolerances." << std::endl;
        return result; // success=false
      }

      // calibration:recalibrate curve: fit the fragment model to all collected peaks, then keep per theoretical fragment
      // (consecutive entries of one m/z and RT) the peak closest to the model, if within the pass window: the errors of
      // the matched peaks, as the other modes measure them
      fragment_curve = MzRecalibration_();
      if (curve)
      {
        if (fragment_errors.size() >= calibration_min_psms_)
        {
          const Size min_peaks = calibration_min_psms_ * fragment_errors.size() / precursor_errors.size();
          // ponytail: one curve per run (a tenth of the peaks does not pin the ends of a cubic: ten curves cost
          // identifications on the validation run); curves per RT chunk with a shared shape if a run's curve drifts
          fragment_curve = MzRecalibration_::fit(fragment_errors, fragment_mzs, fragment_rts, 1, min_peaks, false, true, nullptr, true);
        }
        if (pass == 0 && !fragment_curve.active() && pass_sample > 0)
        {
          sample = 0; // (no curve from the sample: from all spectra of the subset)
          pass = -1;
          continue;
        }
        // (a pass after the first that fits no curve keeps the last one)
        if (pass > 0 && !fragment_curve.active()) { fragment_curve = previous_curve; }
        double moved = std::numeric_limits<double>::infinity();
        if (pass > 0 && fragment_curve.active() && previous_curve.active())
        {
          moved = MzRecalibration_::curveDistance(previous_curve, fragment_curve, fragment_mzs);
          OPENMS_LOG_INFO << "[ProSE] Calibration: pass " << pass + 1 << ", the fragment error curve moved by "
                          << moved << " ppm (95th percentile at the fragments)" << std::endl;
        }
        const bool stable = moved <= curve_tolerance;
        if (fragment_curve.active() && !stable && pass + 1 < max_passes && later_passes)
        {
          // A pass identifies and scores the spectra with the pass window on the errors it sees: the first one on the errors
          // as measured, which holds only the fragments of the best-calibrated m/z range, so fewer and poorer calibration
          // PSMs and a curve that nothing measures again. So the subset is corrected by the last curve, scored again, and
          // the curve is fitted again to the peaks as measured of the PSMs of that pass.
          previous_curve = fragment_curve;
          OPENMS_LOG_INFO << "[ProSE] Calibration: pass " << pass + 2 << " over the subset, corrected by " << fragment_curve.toString() << std::endl;
          correcting = fragment_curve;
          continue;
        }
        else if (pass > 0 && fragment_curve.active())
        {
          OPENMS_LOG_INFO << "[ProSE] Calibration: fragment error curve " << (stable ? "stable" : "not stable") << " after " << pass + 1
                          << " passes" << (stable ? "" : " (the limit)") << std::endl;
        }
        Size kept = 0;
        for (Size i = 0; i < fragment_errors.size();)
        {
          Size j = i, best = i;
          double best_distance = std::numeric_limits<double>::max();
          for (; j < fragment_errors.size() && fragment_mzs[j] == fragment_mzs[i] && fragment_rts[j] == fragment_rts[i]; ++j)
          {
            const double distance = std::fabs(fragment_errors[j] - fragment_curve.error(fragment_mzs[j], fragment_rts[j]));
            if (distance < best_distance) { best_distance = distance; best = j; }
          }
          if (best_distance <= pass_fragment)
          {
            fragment_errors[kept] = fragment_errors[best];
            fragment_mzs[kept] = fragment_mzs[i];
            fragment_rts[kept] = fragment_rts[i];
            ++kept;
          }
          i = j;
        }
        fragment_errors.resize(kept);
        fragment_mzs.resize(kept);
        fragment_rts.resize(kept);
      }
      break;
    }
    if (curve && !fragment_curve.active())
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: no fragment error curve could be fitted (calibration:recalibrate curve); the fragment m/z are not corrected." << std::endl;
    }

    // The share of the errors as measured (before a recalibration) that the measured error recordRecalibration_()
    // reports covers: precursor offset +/- 4 sigma, fragment |offset| + 4 sigma. A measurement that misses more than the
    // windows below accept (two error modes, heavy tails) is recorded as not trusted.
    result.precursor_measured_coverage = measuredCoverage_(precursor_errors, false);
    if (fragment_errors.size() >= calibration_min_psms_) result.fragment_measured_coverage = measuredCoverage_(fragment_errors, true);

    // --- calibration:recalibrate: the error models, which applyRecalibration_() subtracts from the spectra. The
    //     windows below are then derived from the residual errors, the errors the main search sees. ---
    const bool recalibrate = calibration_recalibrate_ != "none";
    if (recalibrate)
    {
      const bool linear = calibration_recalibrate_ == "linear";
      const auto robust = [](std::vector<double> values)
      {
        const auto [med, mad] = medianAndMADBySelection(values);
        return std::make_pair(med, 1.4826 * mad);
      };
      std::tie(result.precursor_raw_shift, result.precursor_raw_sigma) = robust(precursor_errors);
      // (the decoy hits mark the chunks whose PSMs are mostly false matches, which take the model of their neighbours)
      result.precursor_recalibration = MzRecalibration_::fit(precursor_errors, precursor_mzs, precursor_rts,
        calibration_recalibrate_rt_chunks_, calibration_min_psms_, linear, precursor_ppm,
        decoy_string.empty() ? nullptr : &precursor_decoys);
      const auto& chunks = result.precursor_recalibration.chunks;
      const Size from_neighbours = static_cast<Size>(std::count_if(chunks.begin(), chunks.end(),
                                                                   [](const MzRecalibration_::Chunk& c) { return c.from_neighbours; }));
      if (from_neighbours > 0)
      {
        OPENMS_LOG_INFO << "[ProSE] Calibration: in " << from_neighbours << " of " << chunks.size() << " RT chunks of the "
                        << "precursor model most calibration PSMs are false matches (more than a quarter are decoy hits): "
                        << (from_neighbours == chunks.size() ? "one offset for the run, the median of all calibration PSMs."
                                                             : "they take the model of their neighbours.") << std::endl;
      }
      for (Size i = 0; i < precursor_errors.size(); ++i)
      {
        precursor_errors[i] -= result.precursor_recalibration.error(precursor_mzs[i], precursor_rts[i]);
      }
      if (fragment_errors.size() >= calibration_min_psms_)
      {
        std::tie(result.fragment_raw_shift, result.fragment_raw_sigma) = robust(fragment_errors);
        // chunks of as many PSMs as the precursor model's: the peaks of calibration:min_psms PSMs at least
        const Size min_peaks = calibration_min_psms_ * fragment_errors.size() / precursor_errors.size();
        // ('curve' without a ppm tolerance or the peaks as measured: the line)
        result.fragment_recalibration = curve ? fragment_curve : MzRecalibration_::fit(fragment_errors, fragment_mzs, fragment_rts,
          calibration_recalibrate_rt_chunks_, min_peaks, linear || calibration_recalibrate_ == "curve", fragment_mass_tolerance_unit_ppm);
        for (Size i = 0; i < fragment_errors.size(); ++i)
        {
          fragment_errors[i] -= result.fragment_recalibration.error(fragment_mzs[i], fragment_rts[i]);
        }
      }
    }

    // --- The calibrated windows: the median (offset) and 1.4826 * MAD (spread) of the signed errors are robust against
    //     the false matches among the PSMs; quantiles of the errors are not (the accepted false matches spread over the
    //     whole window and make up the tails). Precursor: [median - k sigma, median + k sigma], asymmetric about zero
    //     when the run has an offset; fragment: |median| + k sigma, symmetric (fragment masses are not shifted). Both only
    //     tighten: within [the floor, the configured tolerance]. ---
    const double k = calibration_sigma_factor_;
    const double kp = calibration_precursor_sigma_factor_; // 0: the precursor window is kept
    const double min_tolerance = 1e-6; // avoid non-positive tolerances
    // The share of the calibration PSMs (fragment peaks) a window must cover to be trusted. Precursor errors have heavy
    // tails (5 % of the PSMs of an Astral run beyond 8 ppm at sigma 0.7), and a window that cuts them costs more
    // identifications than it gains (the yield is flat from 5 to 50 ppm): 98 %. Fragment windows gain (10 ppm on
    // Astral: +4 %): 95 %.
    constexpr double min_coverage_precursor = 0.98, min_coverage_fragment = 0.95;
    {
      std::vector<double> errs = precursor_errors; // reordered by the selection
      const auto [med, mad] = medianAndMADBySelection(errs); // = Math::median() / Math::MAD()
      const double sigma = 1.4826 * mad;
      result.precursor_shift = med;
      result.precursor_sigma = sigma;
      // Window of the precursor errors e = observed - theoretical in [-cal_upper, +cal_lower]; with kp 0 the configured
      // window: the calibration then measures and reports the precursor errors and tightens the fragment tolerance only
      // (precursor yields are flat from 5 to 50 ppm; a tightened precursor window cost identifications on every Orbitrap
      // HCD set of the 20-run benchmark at equal entrapment FDP, HF-X -1.6 %, Lumos LFQ -2.8 %)
      const double lo = kp > 0 ? med - kp * sigma : -precursor_mass_tolerance_upper_;
      const double hi = kp > 0 ? med + kp * sigma : precursor_mass_tolerance_lower_;
      const double cal_lower_raw = std::max(min_tolerance, hi);
      const double cal_upper_raw = std::max(min_tolerance, -lo);
      // Zero spread (a synthetic fixture, a single peptide): nothing to calibrate the precursor window on, it is kept;
      // the fragment calibration still applies. Unit-aware: 1 ppm, or min_tolerance in Da.
      result.extreme_bias = kp > 0 && (hi - lo) < (precursor_ppm ? 1.0 : min_tolerance);
      result.precursor_spread = std::max(cal_lower_raw, cal_upper_raw);
      // The median and MAD describe one mode. A run whose offset changes (a lock mass that drops out: half of the
      // precursors near 0 ppm, half near +10 ppm) has errors the window does not cover; then the measurement is not
      // trusted and the configured window stays.
      result.precursor_coverage = std::count_if(precursor_errors.begin(), precursor_errors.end(),
                                                [&](double e) { return e >= lo && e <= hi; }) / static_cast<double>(precursor_errors.size());
      if (kp <= 0)
      {
        result.cal_lower = precursor_mass_tolerance_lower_;
        result.cal_upper = precursor_mass_tolerance_upper_;
      }
      else if (!result.extreme_bias && result.precursor_coverage < min_coverage_precursor)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the precursor window median +/- " << kp << " sigma covers only "
                        << static_cast<int>(100.0 * result.precursor_coverage) << " % of the calibration PSMs (heavy tails, or a mass "
                        << "offset that changes over the run): the configured precursor window is kept." << std::endl;
        result.cal_lower = precursor_mass_tolerance_lower_;
        result.cal_upper = precursor_mass_tolerance_upper_;
      }
      else if (!result.extreme_bias)
      {
        const double floor_side = calibration_min_precursor_tol_ > 0 ? calibration_min_precursor_tol_ : (precursor_ppm ? 2.5 : 0.005);
        result.cal_lower = std::clamp(cal_lower_raw, std::min(floor_side, precursor_mass_tolerance_lower_), precursor_mass_tolerance_lower_);
        result.cal_upper = std::clamp(cal_upper_raw, std::min(floor_side, precursor_mass_tolerance_upper_), precursor_mass_tolerance_upper_);
      }
      // else: cal_lower/cal_upper stay at 0; the writeback skips the precursor window.
    }
    result.psms = precursor_errors.size();
    result.fragment_peaks = fragment_errors.size();
    if (fragment_errors.size() >= calibration_min_psms_)
    {
      const auto [med, mad] = medianAndMADBySelection(fragment_errors);
      const double sigma = 1.4826 * mad;
      result.fragment_shift = med;
      result.fragment_sigma = sigma;
      const double measured = std::fabs(med) + k * sigma;
      // (recalibrated: the residual errors have no offset, so the floor is half)
      const double floor = calibration_min_fragment_tol_ > 0 ? calibration_min_fragment_tol_
                         : (fragment_mass_tolerance_unit_ppm ? (recalibrate ? 5.0 : 10.0) : (recalibrate ? 0.01 : 0.02));
      result.fragment_tolerance = std::clamp(measured, std::min(floor, fragment_mass_tolerance_), fragment_mass_tolerance_);
      result.fragment_coverage = std::count_if(fragment_errors.begin(), fragment_errors.end(),
                                               [&](double e) { return std::fabs(e) <= result.fragment_tolerance; }) / static_cast<double>(fragment_errors.size());
      if (result.fragment_coverage < min_coverage_fragment)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment tolerance |median| + " << k << " sigma covers only "
                        << static_cast<int>(100.0 * result.fragment_coverage) << " % of the matched fragment peaks of the calibration PSMs: "
                        << "the configured fragment tolerance is kept." << std::endl;
        result.fragment_tolerance = fragment_mass_tolerance_;
      }
      if (!fragment_mass_tolerance_unit_ppm && measured < 0.05)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment errors (|offset| + " << k << " sigma = " << std::fixed << std::setprecision(3)
                        << measured << " Da) are those of high-resolution MS2 spectra, but the fragment tolerance is given in Da ("
                        << fragment_mass_tolerance_ << "); a ppm tolerance would suit the data." << std::endl;
      }
    }
    else
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    // Trusted (after the precursor and fragment fits above): the measured error covers the errors (one mode, no heavy
    // tails), and it lies inside the window the calibration pass searched with. A measured range that reaches that window
    // was censored (errors beyond it were not matched, so the median and MAD come from a truncated sample) and is a lower
    // bound, not a tolerance to adopt.
    {
      const double p_off = recalibrate ? result.precursor_raw_shift : result.precursor_shift;
      const double p_sigma = recalibrate ? result.precursor_raw_sigma : result.precursor_sigma;
      const bool p_censored = p_off + 4.0 * p_sigma >= pass_lower || p_off - 4.0 * p_sigma <= -pass_upper;
      result.precursor_measurement_trusted = !p_censored && result.precursor_measured_coverage >= min_coverage_precursor;
      // ('curve': the errors as measured spread over the curve's range and were collected in a wider window; the pass
      // window cuts their distances to the curve, the residual errors)
      const double f_off = recalibrate && !curve ? result.fragment_raw_shift : result.fragment_shift;
      const double f_sigma = recalibrate && !curve ? result.fragment_raw_sigma : result.fragment_sigma;
      const bool f_censored = measuredError_(f_off, f_sigma) >= pass_fragment;
      result.fragment_measurement_trusted = fragment_errors.size() >= calibration_min_psms_ && !f_censored
                                            && result.fragment_measured_coverage >= min_coverage_fragment;
    }

    // A configured window tighter than the error of the run cuts correct matches; the calibration only tightens, so it is
    // reported (a precursor window of 4.5-5 ppm from a methods text or an SDRF is usually a post-recalibration value)
    const double prec_hi = result.precursor_shift + 4.0 * result.precursor_sigma, prec_lo = result.precursor_shift - 4.0 * result.precursor_sigma;
    if (prec_hi > precursor_mass_tolerance_lower_ || -prec_lo > precursor_mass_tolerance_upper_)
    {
      // In the parameters' terms: an error e = observed - theoretical > 0 needs mass_tolerance_lower >= e (the mass-space
      // window is [observed - lower, observed + upper]), e < 0 needs mass_tolerance_upper >= -e
      std::ostringstream sides;
      sides << std::fixed << std::setprecision(3);
      if (prec_hi > precursor_mass_tolerance_lower_)
      {
        sides << " precursor:mass_tolerance_lower " << precursor_mass_tolerance_lower_ << " " << precursor_mass_tolerance_unit_
              << " is below the measured +" << prec_hi << ";";
      }
      if (-prec_lo > precursor_mass_tolerance_upper_)
      {
        sides << " precursor:mass_tolerance_upper " << precursor_mass_tolerance_upper_ << " " << precursor_mass_tolerance_unit_
              << " is below the measured " << -prec_lo << ";";
      }
      OPENMS_LOG_WARN << "[ProSE] Calibration: the precursor window is tighter than the errors the search sees (observed - "
                      << "theoretical from " << std::fixed << std::setprecision(3) << prec_lo << " to " << prec_hi << " "
                      << precursor_mass_tolerance_unit_ << ", offset +/- 4 sigma):" << sides.str() << " this loses correct matches; "
                      << "ProSE gains from wide precursor windows, e.g. 20 ppm" << (recalibrate ? "" : " (or calibration:recalibrate offset)")
                      << "." << std::endl;
    }
    if (fragment_errors.size() >= calibration_min_psms_ && measuredError_(result.fragment_shift, result.fragment_sigma) > fragment_mass_tolerance_)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment tolerance " << fragment_mass_tolerance_ << " " << fragment_mass_tolerance_unit_
                      << " is tighter than the errors the search sees (|offset| + 4 sigma = "
                      << measuredError_(result.fragment_shift, result.fragment_sigma) << ") and loses matched fragments." << std::endl;
    }

    if (!recalibrate)
    {
      result.precursor_raw_shift = result.precursor_shift;
      result.precursor_raw_sigma = result.precursor_sigma;
      result.fragment_raw_shift = result.fragment_shift;
      result.fragment_raw_sigma = result.fragment_sigma;
    }
    result.success = true;

    if (result.recalibrated())
    {
      OPENMS_LOG_INFO << "[ProSE]   Recalibration (" << calibration_recalibrate_ << "): precursor offset " << std::fixed
                      << std::setprecision(2) << result.precursor_raw_shift << ", sigma " << result.precursor_raw_sigma
                      << " -> residual offset " << result.precursor_shift << ", sigma " << result.precursor_sigma << " "
                      << precursor_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.precursor_raw_shift, result.precursor_raw_sigma)
                      << " -> " << measuredError_(result.precursor_shift, result.precursor_sigma) << "); model: " << result.precursor_recalibration.toString() << std::endl;
      if (result.fragment_recalibration.active())
      {
        OPENMS_LOG_INFO << "[ProSE]   Recalibration (" << calibration_recalibrate_ << "): fragment offset " << std::fixed
                        << std::setprecision(2) << result.fragment_raw_shift << ", sigma " << result.fragment_raw_sigma
                        << " -> residual offset " << result.fragment_shift << ", sigma " << result.fragment_sigma << " "
                        << fragment_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.fragment_raw_shift, result.fragment_raw_sigma)
                        << " -> " << measuredError_(result.fragment_shift, result.fragment_sigma) << "); model: " << result.fragment_recalibration.toString() << std::endl;
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Calibration: " << precursor_errors.size() << " PSMs used (top "
                    << static_cast<int>(100.0 * precursor_errors.size() / std::max<Size>(1, cal_hits_total)) << "% by score); pass windows [-"
                    << pass_lower << ", +" << pass_upper << "] " << precursor_mass_tolerance_unit_ << ", fragment "
                    << pass_fragment << " " << fragment_mass_tolerance_unit_ << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor: offset " << std::fixed << std::setprecision(2) << result.precursor_shift
                    << ", sigma " << result.precursor_sigma << " " << precursor_mass_tolerance_unit_
                    << " (|offset| + 4 sigma " << measuredError_(result.precursor_shift, result.precursor_sigma) << ") -> window [-" << result.cal_lower << ", +" << result.cal_upper << "]"
                    << " (covers " << static_cast<int>(100.0 * result.precursor_coverage) << " % of the PSMs)"
                    << (kp <= 0 ? " (the configured window, calibration:precursor_sigma_factor 0)" : "")
                    << (result.extreme_bias ? " (zero spread: the configured window is kept)" : "") << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Fragment:  offset " << result.fragment_shift << ", sigma " << result.fragment_sigma
                    << " " << fragment_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.fragment_shift, result.fragment_sigma)
                    << "; " << result.fragment_peaks << " peaks) -> tolerance "
                    << fragment_mass_tolerance_ << " -> " << result.fragment_tolerance << " " << fragment_mass_tolerance_unit_
                    << " (covers " << static_cast<int>(100.0 * result.fragment_coverage) << " % of the peaks)" << std::endl;

    return result;
  }

  // =====================================================================
  // calibration:recalibrate: the m/z error models and their application to the spectra
  // =====================================================================
  ProSEAlgorithm::MzRecalibration_::AtRT ProSEAlgorithm::MzRecalibration_::at(double rt) const
  {
    AtRT model;
    model.ppm = ppm;
    if (chunks.empty()) return model;
    if (chunks.size() == 1 || rt <= chunks.front().rt) { model.lo = &chunks.front(); return model; }
    if (rt >= chunks.back().rt) { model.lo = &chunks.back(); return model; }
    // the chunks around rt: lo->rt <= rt < hi->rt
    const auto hi = std::upper_bound(chunks.begin(), chunks.end(), rt, [](double r, const Chunk& c) { return r < c.rt; });
    const auto lo = hi - 1;
    model.lo = &*lo;
    model.hi = &*hi;
    model.w = (rt - lo->rt) / (hi->rt - lo->rt);
    return model;
  }

  double ProSEAlgorithm::MzRecalibration_::error(double mz, double rt) const
  {
    return at(rt).error(mz);
  }

  double ProSEAlgorithm::MzRecalibration_::correct(double mz, double rt) const
  {
    return at(rt).correct(mz);
  }

  double ProSEAlgorithm::MzRecalibration_::uncorrect(double corrected, double rt) const
  {
    if (chunks.empty()) return corrected;
    // correct() maps mz to mz (1 - e(mz) 1e-6), or mz - e(mz) in Da: iterate mz = corrected / (1 - e(mz) 1e-6). Exact
    // after one step for 'offset' (its error does not depend on m/z); for a line in m/z each step shrinks the
    // remaining difference by a factor of about slope * mz * 1e-6 (a curve: its slope), so three steps suffice.
    double mz = corrected;
    for (int step = 0; step < 3; ++step)
    {
      const double e = error(mz, rt);
      const double next = ppm ? corrected / (1.0 - e * 1e-6) : corrected + e;
      if (next == mz) break; // converged
      mz = next;
    }
    return mz;
  }

  double ProSEAlgorithm::MzRecalibration_::curveDistance(const MzRecalibration_& a, const MzRecalibration_& b,
                                                         const std::vector<double>& mzs)
  {
    if (!a.curve() || !b.curve() || mzs.empty()) return std::numeric_limits<double>::infinity();
    const Chunk& ca = a.chunks.front();
    const Chunk& cb = b.chunks.front();
    std::vector<double> differences(mzs.size());
    for (Size k = 0; k < mzs.size(); ++k) { differences[k] = std::fabs(ca.error(mzs[k]) - cb.error(mzs[k])); }
    const auto q95 = differences.begin() + static_cast<std::ptrdiff_t>(0.95 * static_cast<double>(differences.size() - 1));
    std::nth_element(differences.begin(), q95, differences.end());
    return *q95;
  }

  std::string ProSEAlgorithm::MzRecalibration_::toString() const
  {
    if (chunks.empty()) return "";
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << (ppm ? "ppm" : "Da");
    for (const Chunk& c : chunks)
    {
      s << "; rt " << std::fixed << std::setprecision(1) << c.rt << ": ";
      if (c.s_hi > 0.0)
      {
        s << "cubic in t = sqrt(mz) scaled from [" << std::setprecision(4) << c.s_lo << ", " << c.s_hi << "] to [-1, 1]: "
          << c.intercept << " " << c.slope << " " << c.c2 << " " << c.c3 << " (" << c.points << ")";
        continue;
      }
      s << std::setprecision(4) << c.intercept;
      if (c.slope != 0.0) { s << (c.slope < 0 ? " - " : " + ") << std::scientific << std::setprecision(3) << std::fabs(c.slope) << " mz"; }
      s << " (" << c.points;
      if (c.from_neighbours) { s << ", " << c.decoys << " decoys: mostly false matches"; }
      s << ")";
    }
    return s.str();
  }

  ProSEAlgorithm::MzRecalibration_ ProSEAlgorithm::MzRecalibration_::fit(const std::vector<double>& errors, const std::vector<double>& mzs,
    const std::vector<double>& rts, Size n_chunks, Size min_points, bool linear, bool ppm, const std::vector<char>* decoys, bool curve)
  {
    MzRecalibration_ model;
    model.ppm = ppm;
    const Size n = errors.size();
    if (n == 0 || mzs.size() != n || rts.size() != n || (decoys != nullptr && decoys->size() != n)) return model;
    for (Size i = 0; i < n; ++i) // (values that are not numbers break the orderings below and the histograms of fitErrorCurve())
    {
      if (!std::isfinite(errors[i]) || !std::isfinite(mzs[i]) || !std::isfinite(rts[i])) return model;
    }
    std::vector<Size> order(n);
    std::iota(order.begin(), order.end(), Size(0));
    std::stable_sort(order.begin(), order.end(), [&rts](Size a, Size b) { return rts[a] < rts[b]; });
    n_chunks = std::max<Size>(1, std::min(n_chunks, n / std::max<Size>(1, min_points)));
    std::vector<char> mostly_false(n_chunks, 0); // (with decoys only; never with curve, whose chunk may be left out)
    for (Size c = 0; c < n_chunks; ++c)
    {
      const Size begin = n * c / n_chunks, end = n * (c + 1) / n_chunks;
      std::vector<double> e, x;
      e.reserve(end - begin);
      x.reserve(end - begin);
      for (Size i = begin; i < end; ++i) { e.push_back(errors[order[i]]); x.push_back(mzs[order[i]]); }
      std::vector<double> sorted = e;
      const auto [med, mad] = medianAndMADBySelection(sorted);
      Chunk chunk{rts[order[(begin + end) / 2]], med, 0.0, e.size()}; // the centre: the median RT of the chunk
      if (decoys != nullptr)
      {
        for (Size i = begin; i < end; ++i) { chunk.decoys += (*decoys)[order[i]] != 0; }
        mostly_false[c] = 4 * chunk.decoys > end - begin; // estimated false matches (2 * decoys) > half of the points
      }
      if (curve && ppm)
      {
        Cubic f;
        Size points = 0;
        if (fitErrorCurve(e, x, f, points))
        {
          chunk.intercept = f.c[0];
          chunk.slope = f.c[1];
          chunk.c2 = f.c[2];
          chunk.c3 = f.c[3];
          chunk.s_lo = f.lo;
          chunk.s_hi = f.hi;
          chunk.points = points;
        }
        else { continue; } // no chunk: the median of all peaks within the collection window is no calibration
      }
      else if (linear || curve)
      {
        // the least-squares line through the points within 3 sigma of the median (the false matches among the PSMs)
        const double cut = 3.0 * 1.4826 * mad;
        std::vector<double> te, tx;
        for (Size i = 0; i < e.size(); ++i) { if (std::fabs(e[i] - med) <= cut) { te.push_back(e[i]); tx.push_back(x[i]); } }
        MZTrafoModel line(ppm);
        double intercept = 0.0, slope = 0.0, power = 0.0;
        if (te.size() >= 2 && line.train(te, tx, std::vector<double>(te.size(), 1.0), MZTrafoModel::MODELTYPE::LINEAR, false))
        {
          line.getCoefficients(intercept, slope, power);
          if (std::isfinite(intercept) && std::isfinite(slope)) { chunk.intercept = intercept; chunk.slope = slope; chunk.points = te.size(); }
        }
      }
      model.chunks.push_back(chunk);
    }
    // A chunk whose calibration PSMs are mostly false matches takes the model of its neighbours. Its median is that of
    // the false matches, which need not be near the mass error of the instrument: a contaminant that elutes over
    // minutes (PEG late in a run) is matched to targets and decoys alike at the same +11 ppm, and the last three of
    // eight chunks of a PXD019643 class I run (twice their decoys: 78-96 % of their PSMs) got +9 to +11 ppm, whose
    // correct PSMs are within 0.3 ppm of 0, and lost them from the 5 ppm window of the search (OpenMS/OpenMS#10885).
    // The decoys tell such a chunk from a drift of the instrument, which shifts the errors of correct PSMs and is
    // followed (-15 ppm in the last chunk of a PXD011695 EThcD run, 1 decoy among its 88 PSMs).
    if (std::find(mostly_false.begin(), mostly_false.end(), 1) == mostly_false.end()) return model;
    if (std::find(mostly_false.begin(), mostly_false.end(), 0) == mostly_false.end())
    {
      // no chunk to take the model from: the median of all points, one offset for the run
      std::vector<double> all = errors;
      const double median = medianAndMADBySelection(all).first;
      for (Chunk& chunk : model.chunks) { chunk.intercept = median; chunk.slope = 0.0; chunk.from_neighbours = true; }
      return model;
    }
    const std::vector<Chunk> fitted = model.chunks;
    for (Size c = 0; c < n_chunks; ++c)
    {
      if (!mostly_false[c]) continue;
      SignedSize before = static_cast<SignedSize>(c) - 1;
      while (before >= 0 && mostly_false[static_cast<Size>(before)]) { --before; }
      Size after = c + 1;
      while (after < n_chunks && mostly_false[after]) { ++after; }
      Chunk& chunk = model.chunks[c];
      chunk.from_neighbours = true;
      if (before < 0 || after == n_chunks || !(fitted[after].rt > fitted[static_cast<Size>(before)].rt))
      {
        const Chunk& nearest = before < 0 ? fitted[after] : fitted[static_cast<Size>(before)]; // at an end (or no RT range): the nearest
        chunk.intercept = nearest.intercept;
        chunk.slope = nearest.slope;
        continue;
      }
      const Chunk& lo = fitted[static_cast<Size>(before)];
      const Chunk& hi = fitted[after];
      const double w = (chunk.rt - lo.rt) / (hi.rt - lo.rt); // as the model interpolates between chunks
      chunk.intercept = (1.0 - w) * lo.intercept + w * hi.intercept;
      chunk.slope = (1.0 - w) * lo.slope + w * hi.slope;
    }
    return model;
  }

  void ProSEAlgorithm::applyRecalibration_(const CalibrationResult_& cal, PeakMap& spectra, PeakMap* query_spectra,
                                           PeakMap* evidence_spectra, FragmentIonLikelihoodModel::PeakLists* ion_evidence,
                                           PeakCharges_* charges) const
  {
    if (!cal.success || !cal.recalibrated()) return;
    const MzRecalibration_& precursor_model = cal.precursor_recalibration;
    const MzRecalibration_& fragment_model = cal.fragment_recalibration;
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt observed_key = registry.registerName(OBSERVED_PRECURSOR_MZ);
    // the co-isolated precursors of a chimeric search (CoIsolationDetector::annotate(); UInt(-1): none annotated in this
    // process), from the same survey scans as the selected precursor: corrected by its model
    const UInt coisolated_key = registry.getIndex(CoIsolationDetector::META_MZ);
    const UInt observed_coisolated_key = coisolated_key != UInt(-1) ? registry.registerName(OBSERVED_COISOLATED_MZ) : UInt(-1);
    // The model at the spectrum's RT once, per peak only its arithmetic (the chunk lookup per peak was a third of the
    // time of this function in stack samples of a 27-run search)
    const auto correct_peaks = [](MSSpectrum& s, const MzRecalibration_::AtRT& model)
    {
      for (Peak1D& p : s) { p.setMZ(model.correct(p.getMZ())); } // (monotone for any sane model: the order stays)
    };
    // 'curve' gives the error at the m/z a peak was measured at, which a deisotoped peak no longer has (its m/z is that
    // of the singly charged ion): each peak is corrected at its measured m/z (in place, below), or, without the charges
    // of the peaks, the searched spectra are made again from the corrected peaks as measured (fragment:query_spectrum
    // raw; a curve is only fitted with these peaks at hand).
    const bool curve = fragment_model.curve() && query_spectra != nullptr;
    // In place with the charges of the deisotoped peaks: each peak at the m/z it was measured at, (m + (z - 1) proton) / z,
    // then moved back to its singly charged ion (as Deisotoper::deisotopeAndSingleCharge()). The ion priors of all peaks
    // are taken from the corrected fragment-evidence copies, which hold the same peaks; without these copies, rebuilt.
    // (the ion priors' lists exist with annotate:self_trained_ion_priors only, as preprocessSpectra_() fills them)
    FragmentIonLikelihoodModel::PeakLists* const ion_lists
      = self_trained_ion_priors_ && ion_evidence != nullptr && ion_evidence->size() == spectra.size() ? ion_evidence : nullptr;
    const bool in_place = curve && charges != nullptr && charges->searched.size() == spectra.size()
                          && (evidence_spectra == nullptr || charges->deisotoped.size() == spectra.size())
                          && (ion_lists == nullptr || ion_prior_scored_peaks_ || evidence_spectra != nullptr);
    const bool rebuild = curve && !in_place;
    Size precursors = 0, coisolated = 0, peaks = 0;
#pragma omp parallel for schedule(static) reduction(+: precursors, coisolated, peaks)
    for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
    {
      MSSpectrum& s = spectra[i];
      if (s.getMSLevel() != 2) continue;
      const double rt = s.getRT();
      if (precursor_model.active() && !s.getPrecursors().empty())
      {
        Precursor& p = s.getPrecursors()[0];
        if (!p.metaValueExists(observed_key)) { p.setMetaValue(observed_key, p.getMZ()); } // (a second correction keeps the first observation)
        const double corrected = precursor_model.correct(p.getMZ(), rt);
        p.setMZ(corrected);
        if (query_spectra != nullptr && !(*query_spectra)[i].getPrecursors().empty()) { (*query_spectra)[i].getPrecursors()[0].setMZ(corrected); }
        ++precursors;
      }
      if (precursor_model.active() && coisolated_key != UInt(-1) && s.metaValueExists(coisolated_key))
      {
        DoubleList mzs = s.getMetaValue(coisolated_key);
        if (!s.metaValueExists(observed_coisolated_key)) { s.setMetaValue(observed_coisolated_key, mzs); } // (as above)
        for (double& mz : mzs) { mz = precursor_model.correct(mz, rt); }
        coisolated += mzs.size();
        s.setMetaValue(coisolated_key, std::move(mzs));
      }
      if (fragment_model.active())
      {
        const MzRecalibration_::AtRT model = fragment_model.at(rt);
        if (query_spectra != nullptr) { correct_peaks((*query_spectra)[i], model); }
        if (in_place)
        {
          correctAtMeasuredMZ(s, model, charges->searched[i]);
          if (evidence_spectra != nullptr) { correctAtMeasuredMZ((*evidence_spectra)[i], model, charges->deisotoped[i]); }
          if (ion_lists != nullptr) { ion_lists->assign(i, ion_prior_scored_peaks_ ? s : (*evidence_spectra)[i]); }
        }
        else if (rebuild)
        {
          const MSSpectrum& measured = (*query_spectra)[i];
          // ponytail: clear() drops the data arrays too; those preprocessSpectra_() makes are made again, arrays of the
          // input carried along (e.g. ion mobility) are lost: copy them from the measured spectrum if a later step reads them
          s.clear(false); // the peaks; the settings (RT, precursors) stay
          s.insert(s.end(), measured.begin(), measured.end());
        }
        else
        {
          correct_peaks(s, model);
          if (evidence_spectra != nullptr) { correct_peaks((*evidence_spectra)[i], model); }
        }
        peaks += s.size();
      }
    }
    if (rebuild) // (the copies are made again with the spectra)
    {
      preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm", deisotope_requested_,
                         peaks_keep_n_, peaks_window_top_, peaks_window_type_, deisotoping_,
                         self_trained_ion_priors_ ? ion_evidence : nullptr, ion_prior_scored_peaks_, evidence_spectra,
                         query_spectra);
    }
    if (charges != nullptr) { *charges = PeakCharges_(); } // (only needed here)
    OPENMS_LOG_INFO << "[ProSE] Recalibration (" << calibration_recalibrate_ << "): corrected the precursor m/z of " << precursors
                    << " spectra" << (coisolated > 0 ? " and of " + std::to_string(coisolated) + " co-isolated precursors" : std::string())
                    << " and " << peaks << " fragment peaks in place (the PSMs report the observed m/z)" << std::endl;
  }

  double ProSEAlgorithm::measuredCoverage_(const std::vector<double>& errors, bool fragment)
  {
    if (errors.empty()) return 0.0;
    std::vector<double> values = errors; // reordered by the selection
    const std::pair<double, double> med_mad = medianAndMADBySelection(values);
    const double med = med_mad.first; // (no structured binding: clang cannot capture one in a lambda with OpenMP)
    const double sigma = 1.4826 * med_mad.second;
    const double measured = measuredError_(med, sigma);
    const auto covered = [&](double e) { return fragment ? std::fabs(e) <= measured : (e >= med - 4.0 * sigma && e <= med + 4.0 * sigma); };
    return std::count_if(errors.begin(), errors.end(), covered) / static_cast<double>(errors.size());
  }

  void ProSEAlgorithm::recordRecalibration_(const CalibrationResult_& cal, std::vector<ProteinIdentification>& protein_ids) const
  {
    if (protein_ids.empty() || !cal.success) return;
    ProteinIdentification::SearchParameters sp = protein_ids[0].getSearchParameters();
    // The mass errors the calibration pass measured: offset (median), sigma (1.4826 * MAD) and the measured error
    // |offset| + 4 sigma, the range that holds practically all correct matches (in the tolerance units)
    sp.setMetaValue("calibration:precursor_error_offset", cal.precursor_raw_shift);
    sp.setMetaValue("calibration:precursor_error_sigma", cal.precursor_raw_sigma);
    sp.setMetaValue("calibration:precursor_measured_error", measuredError_(cal.precursor_raw_shift, cal.precursor_raw_sigma));
    // (no fragment records without a fragment measurement, fewer than calibration:min_psms matched peaks: their zeros
    // would read as a perfect fragment calibration)
    const bool fragments_measured = cal.fragment_peaks >= calibration_min_psms_;
    if (fragments_measured)
    {
      sp.setMetaValue("calibration:fragment_error_offset", cal.fragment_raw_shift);
      sp.setMetaValue("calibration:fragment_error_sigma", cal.fragment_raw_sigma);
      sp.setMetaValue("calibration:fragment_measured_error", measuredError_(cal.fragment_raw_shift, cal.fragment_raw_sigma));
    }
    // Whether a measurement describes the run: it covers at least 98 % (precursor) or 95 % (fragment peaks) of the
    // calibration PSMs, lies inside the calibration pass's window (not censored), and the fragments come from
    // calibration:min_psms peaks or more. A pipeline that takes the measured error as a tolerance should not take an
    // untrusted one.
    sp.setMetaValue("calibration:psms", static_cast<Int>(cal.psms));
    sp.setMetaValue("calibration:precursor_measurement_trusted", cal.precursor_measurement_trusted ? "true" : "false");
    sp.setMetaValue("calibration:fragment_measurement_trusted", cal.fragment_measurement_trusted ? "true" : "false");
    if (cal.recalibrated())
    {
      sp.setMetaValue("calibration:recalibrate", calibration_recalibrate_);
      sp.setMetaValue("calibration:recalibrate_rt_chunks", static_cast<Int>(calibration_recalibrate_rt_chunks_));
      sp.setMetaValue("calibration:precursor_model", cal.precursor_recalibration.toString());
      sp.setMetaValue("calibration:precursor_residual_offset", cal.precursor_shift);
      sp.setMetaValue("calibration:precursor_residual_sigma", cal.precursor_sigma);
      sp.setMetaValue("calibration:precursor_residual_measured_error", measuredError_(cal.precursor_shift, cal.precursor_sigma));
      if (fragments_measured)
      {
        sp.setMetaValue("calibration:fragment_model", cal.fragment_recalibration.toString());
        sp.setMetaValue("calibration:fragment_residual_offset", cal.fragment_shift);
        sp.setMetaValue("calibration:fragment_residual_sigma", cal.fragment_sigma);
        sp.setMetaValue("calibration:fragment_residual_measured_error", measuredError_(cal.fragment_shift, cal.fragment_sigma));
      }
    }
    protein_ids[0].setSearchParameters(std::move(sp));
  }

  // =====================================================================
  // Serialize the end-of-search report to a YAML string. Hand-rolled (no YAML
  // library), so neither the TOPP tool nor the library needs an extra dependency
  // for this. Every string scalar is double-quoted (q()) so values containing ':'
  // (e.g. Windows paths) or other metacharacters cannot break the structure, and
  // non-finite numbers are emitted as null (num()).
  // =====================================================================
  std::string ProSEAlgorithm::renderRunSummaryYaml(
      const MultiFileSearchResult& mfres,
      const std::vector<std::pair<std::string, std::vector<std::string>>>& manifest,
      Size files_failed,
      Size files_total)
  {
    // Double-quoted YAML scalar: escape the two structural chars, the common
    // whitespace escapes, and every other C0 control + DEL as \xNN. A literal
    // control character is ill-formed in a double-quoted scalar and would corrupt
    // or invalidate the document — reachable via arbitrary input file paths.
    auto q = [](const std::string& s) {
      static const char* const hex = "0123456789ABCDEF";
      std::string out;
      out.reserve(s.size() + 2);
      out += '"';
      for (char ch : s)
      {
        const unsigned char c = static_cast<unsigned char>(ch);
        if (c == '\\' || c == '"') { out += '\\'; out += static_cast<char>(c); }
        else if (c == '\n') { out += "\\n"; }
        else if (c == '\t') { out += "\\t"; }
        else if (c == '\r') { out += "\\r"; }
        else if (c < 0x20 || c == 0x7F) { out += "\\x"; out += hex[c >> 4]; out += hex[c & 0x0F]; }
        else { out += static_cast<char>(c); }
      }
      out += '"';
      return out;
    };
    // Locale-independent number; non-finite -> null. YAML 1.1 reads an exponent
    // without a mantissa dot ("8e-05") as a string, so ensure a '.' before 'e'.
    auto num = [](double v) -> std::string {
      if (!std::isfinite(v)) { return "null"; }
      std::ostringstream o;
      o.imbue(std::locale::classic());
      o << v;
      std::string s = o.str();
      const auto e = s.find_first_of("eE");
      if (e != std::string::npos && s.find('.') == std::string::npos) { s.insert(e, ".0"); }
      return s;
    };
    auto bol = [](bool b) -> std::string { return b ? "true" : "false"; };

    std::ostringstream y;
    const SharedSearchStats& sh = mfres.shared;

    y << "shared:\n";
    y << "  database_file: " << q(sh.database_file) << "\n";
    y << "  enzyme: " << q(sh.enzyme) << "\n";
    y << "  precursor_tol_lower: " << num(sh.precursor_tol_lower) << "\n";
    y << "  precursor_tol_upper: " << num(sh.precursor_tol_upper) << "\n";
    y << "  precursor_tol_unit: " << q(sh.precursor_tol_unit) << "\n";
    y << "  fragment_tol: " << num(sh.fragment_tol) << "\n";
    y << "  fragment_tol_unit: " << q(sh.fragment_tol_unit) << "\n";
    y << "  min_charge: " << sh.min_charge << "\n";
    y << "  max_charge: " << sh.max_charge << "\n";
    y << "  missed_cleavages: " << sh.missed_cleavages << "\n";
    auto str_list = [&](const char* key, const std::vector<std::string>& items) {
      y << "  " << key << ":";
      if (items.empty()) { y << " []\n"; return; }
      y << "\n";
      for (const auto& it : items) { y << "    - " << q(it) << "\n"; }
    };
    str_list("fixed_mods", std::vector<std::string>(sh.fixed_mods.begin(), sh.fixed_mods.end()));
    str_list("variable_mods", std::vector<std::string>(sh.variable_mods.begin(), sh.variable_mods.end()));
    str_list("ion_series", std::vector<std::string>(sh.ion_series.begin(), sh.ion_series.end()));
    y << "  open_search: " << bol(sh.open_search) << "\n";
    y << "  calibration_enabled: " << bol(sh.calibration_enabled) << "\n";
    y << "  peptide_mass_mode: " << bol(sh.peptide_mass_mode) << "\n";
    y << "  chunked: " << bol(sh.chunked) << "\n";
    y << "  decoy_mode: " << q(sh.decoy_mode) << "\n";
    y << "  psm_fdr_threshold: " << num(sh.psm_fdr_threshold) << "\n";
    y << "  protein_fdr_threshold: " << num(sh.protein_fdr_threshold) << "\n";
    y << "  db_target_proteins: " << sh.db_target_proteins << "\n";
    y << "  db_decoy_proteins: " << sh.db_decoy_proteins << "\n";
    y << "  indexed_peptides: " << sh.indexed_peptides << "\n";
    y << "  indexed_fragments: " << sh.indexed_fragments << "\n";
    y << "  seconds_index_build: " << num(sh.seconds_index_build) << "\n";
    y << "  seconds_total: " << num(sh.seconds_total) << "\n";
    y << "  seconds_output: " << num(sh.seconds_output) << "\n";
    y << "  seconds_qpx_extended: " << num(sh.seconds_qpx_extended) << "\n";

    y << "per_file:";
    if (mfres.per_file.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& pf : mfres.per_file)
      {
        const RunStatistics& st = pf.stats;
        // The first key of each list item carries the "- " marker; the rest indent to 4.
        y << "  - input_file: " << q(st.input_file) << "\n";
        y << "    exit_code: " << static_cast<int>(pf.exit_code) << "\n";
        y << "    ms2_spectra: " << st.ms2_spectra << "\n";
        y << "    profile_spectra: " << st.profile_spectra << "\n";
        y << "    charge_out_of_range_spectra: " << st.charge_out_of_range_spectra << "\n";
        y << "    matched_spectra: " << st.matched_spectra << "\n";
        y << "    target_psms: " << st.target_psms << "\n";
        y << "    decoy_psms: " << st.decoy_psms << "\n";
        y << "    fdr_applied: " << bol(st.fdr_applied) << "\n";
        y << "    achieved_psm_fdr: " << num(st.achieved_psm_fdr) << "\n";
        y << "    unique_peptides: " << st.unique_peptides << "\n";
        y << "    unique_proteins: " << st.unique_proteins << "\n";
        y << "    hyperscore:\n";
        y << "      valid: " << bol(st.score_stats_valid) << "\n";
        y << "      min: " << num(st.hyperscore_min) << "\n";
        y << "      median: " << num(st.hyperscore_median) << "\n";
        y << "      max: " << num(st.hyperscore_max) << "\n";
        y << "    charge_histogram:";
        if (st.charge_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [z, c] : st.charge_histogram) { y << "      " << q(std::to_string(z)) << ": " << c << "\n"; } }
        y << "    missed_cleavage_histogram:";
        if (st.missed_cleavage_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [m, c] : st.missed_cleavage_histogram) { y << "      " << q(std::to_string(m)) << ": " << c << "\n"; } }
        y << "    precursor_error:\n";
        y << "      valid: " << bol(st.prec_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.prec_err_median) << "\n";
        y << "      mad_ppm: " << num(st.prec_err_mad) << "\n";
        y << "      recommended_ppm: " << num(st.prec_err_recommended) << "\n";
        y << "    fragment_error:\n";
        y << "      valid: " << bol(st.frag_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.frag_err_median) << "\n";
        y << "      mad_ppm: " << num(st.frag_err_mad) << "\n";
        if (st.frag_err_after_curve) { y << "      after_curve: true\n"; } // (only then: the other reports stay as they were)
        y << "      recommended_ppm: " << (st.frag_err_after_curve ? std::string("null") : num(st.frag_err_recommended)) << "\n";
        y << "    timing_seconds:\n";
        y << "      calibration: " << num(st.seconds_calibration) << "\n";
        y << "      search: " << num(st.seconds_search) << "\n";
        y << "      fdr: " << num(st.seconds_fdr) << "\n";
        if (!st.recalibration.empty()) // calibration: the measured errors; with recalibrate the mode, models and residuals
        {
          y << "    calibration:\n";
          for (const auto& [name, value] : st.recalibration)
          {
            y << "      " << name << ": ";
            if (value.valueType() == DataValue::DOUBLE_VALUE) { y << num(static_cast<double>(value)); }
            else if (value.valueType() == DataValue::INT_VALUE) { y << static_cast<long>(value); }
            else { y << q(value.toString()); }
            y << "\n";
          }
        }
      }
    }

    y << "outputs:";
    if (manifest.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& [label, paths] : manifest)
      {
        y << "  - type: " << q(label) << "\n";
        y << "    paths:";
        if (paths.empty()) { y << " []\n"; }
        else { y << "\n"; for (const auto& p : paths) { y << "      - " << q(p) << "\n"; } }
      }
    }
    y << "files_failed: " << files_failed << "\n";
    y << "files_total: " << files_total << "\n";

    return y.str();
  }

  // =====================================================================
  // Render the modification-discovery section (open search) to an ostream.
  // =====================================================================
  void ProSEAlgorithm::renderModificationSummary(
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      std::ostream& os)
  {
    const auto& dm_stats = mod_analysis.delta_mass_stats;
    const auto& ptm_stats = mod_analysis.ptm_stats;

    os << "[ProSE] ---------------- Modification discovery ---------------\n";
    os << "[ProSE] Delta mass : " << dm_stats.modified_psms << " modified / "
       << dm_stats.total_psms << " PSMs ("
       << std::fixed << std::setprecision(1)
       << (dm_stats.total_psms > 0 ? (100.0 * dm_stats.modified_psms / dm_stats.total_psms) : 0.0)
       << "%), median Δ=" << std::setprecision(4) << dm_stats.median_delta_mass
       << " Da, " << dm_stats.entries.size() << " bins\n";
    os << "[ProSE] PTMs      : " << ptm_stats.total_modified_psms << " PSMs with known PTMs, "
       << ptm_stats.unknown_modification_psms << " unknown, "
       << ptm_stats.num_unique_modifications << " unique PTMs\n";

    if (!ptm_stats.entries.empty())
    {
      os << "[ProSE]   Top PTMs (rank | name | count | % | mass Da):\n";
      size_t rank = 1;
      for (const auto& ptm : ptm_stats.entries)
      {
        if (rank > 15) { break; }
        std::string name = ptm.name;
        if (name.size() > 30) { name = StringUtils::substr(name, 0, 27) + "..."; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | "
           << std::setw(31) << std::left << name << std::right << " | "
           << std::setw(6) << ptm.count << " | "
           << std::setw(5) << std::fixed << std::setprecision(1) << ptm.percentage << " | "
           << std::setw(9) << std::fixed << std::setprecision(4) << ptm.theoretical_mass << "\n";
      }
    }

    std::vector<OpenSearchModificationAnalysis::DeltaMassEntry> unknown_dm;
    for (const auto& entry : dm_stats.entries)
    {
      if (!entry.is_known_modification && entry.count >= 5) { unknown_dm.push_back(entry); }
    }
    if (!unknown_dm.empty())
    {
      os << "[ProSE]   Top unknown delta masses (potential novel PTMs):\n";
      size_t rank = 1;
      for (const auto& dm : unknown_dm)
      {
        if (rank > 10) { break; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | Δ="
           << std::setw(11) << std::fixed << std::setprecision(4) << dm.delta_mass << " Da | "
           << std::setw(6) << dm.count << " PSMs | "
           << dm.unique_peptides << " peptides\n";
      }
    }
  }

  // =====================================================================
  // Render a human-readable single-run summary block.
  // =====================================================================
  std::string ProSEAlgorithm::renderTotalTime(const SharedSearchStats& sh)
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << "total " << sh.seconds_total << "s | search of the file(s) "
       << std::max(0.0, sh.seconds_total - sh.seconds_output) << "s | outputs " << sh.seconds_output << "s";
    if (sh.seconds_qpx_extended > 0.0) { os << " (of which -qpx_extended " << sh.seconds_qpx_extended << "s)"; }
    return os.str();
  }

  void ProSEAlgorithm::renderRunSummary(
      const RunStatistics& s,
      const SharedSearchStats& sh,
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      bool is_open_search,
      std::ostream& os)
  {
    auto join = [](const std::vector<std::string>& v) -> std::string {
      if (v.empty()) { return "(none)"; }
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? ", " : "") + v[i]; }
      return out;
    };
    auto join_int = [](const std::vector<std::string>& v) -> std::string {
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? "," : "") + v[i]; }
      return out.empty() ? std::string("(none)") : out;
    };

    if (!s.input_file.empty())
    {
      os << "[ProSE] Input        : " << s.input_file << "  (" << s.ms2_spectra << " MS2 spectra)\n";
    }
    if (s.ms2_spectra == 0)
    {
      os << "[ProSE] Warning      : no MS2 spectra in the input (e.g. an MS1-only file): nothing was searched\n";
    }
    if (s.profile_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.profile_spectra << " MS2 spectra are profile spectra; centroid them first\n";
    }
    if (s.charge_out_of_range_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.charge_out_of_range_spectra << " MS2 spectra not searched: precursor charge outside of "
         << sh.min_charge << "-" << sh.max_charge << "\n";
    }
    // -- Configuration recap --
    os << "[ProSE] Config       : " << sh.enzyme << ", " << sh.missed_cleavages << " MC | prec [-"
       << sh.precursor_tol_lower << ", +" << sh.precursor_tol_upper << "] " << sh.precursor_tol_unit
       << " | frag " << sh.fragment_tol << " " << sh.fragment_tol_unit
       << " | z " << sh.min_charge << "-" << sh.max_charge << "\n";
    os << "[ProSE]                fixed: " << join(sh.fixed_mods)
       << " | variable: " << join(sh.variable_mods) << "\n";
    os << "[ProSE]                ions: " << join_int(sh.ion_series)
       << " | calibration: " << (sh.calibration_enabled ? "on" : "off")
       << (sh.recalibrate.empty() || sh.recalibrate == "none" ? std::string() : " (recalibrate " + sh.recalibrate + ")")
       << " | mode: " << (sh.open_search ? "open" : "closed")
       << (sh.peptide_mass_mode ? " | low_memory" : "")
       << (sh.chunked ? " | chunked" : "") << "\n";

    // -- Database / fragment index --
    os << "[ProSE] Database     : " << (sh.database_file.empty() ? std::string("(in-memory)") : std::string(sh.database_file))
       << "  (" << sh.db_target_proteins << " target";
    if (sh.db_decoy_proteins > 0) { os << " + " << sh.db_decoy_proteins << " decoy"; }
    os << " proteins, decoys: " << (sh.decoy_mode.empty() ? std::string("n/a") : std::string(sh.decoy_mode)) << ")\n";
    os << "[ProSE] Fragment idx : " << sh.indexed_peptides << " peptides / "
       << sh.indexed_fragments << " fragments  (build "
       << std::fixed << std::setprecision(1) << sh.seconds_index_build << " s)\n";

    // -- Results --
    os << "[ProSE] --------------------------- Results ---------------------------\n";
    os << "[ProSE] Matched      : " << s.matched_spectra << " / " << s.ms2_spectra << " spectra";
    if (s.ms2_spectra > 0)
    {
      os << "   (ID rate " << std::fixed << std::setprecision(1)
         << (100.0 * s.matched_spectra / s.ms2_spectra) << "%)";
    }
    os << "\n";
    os << "[ProSE] Peptides/prot: " << s.unique_peptides << " unique peptides | "
       << s.unique_proteins << " unique proteins\n";

    // FDR / target-decoy. Decide target-only strictly from the database, not from
    // a zero PSM count (a decoy DB with no decoy PSMs is still FDR-capable).
    if (sh.decoy_mode == "none (target-only)")
    {
      os << "[ProSE] FDR          : n/a (target-only database)\n";
    }
    else
    {
      os << "[ProSE] Target/decoy : " << s.target_psms << " target / " << s.decoy_psms << " decoy PSMs";
      if (s.fdr_applied)
      {
        os << "  | PSM FDR <= " << std::fixed << std::setprecision(1) << (sh.psm_fdr_threshold * 100.0) << "%";
        if (s.achieved_psm_fdr >= 0.0)
        {
          os << ", achieved " << std::setprecision(2) << (s.achieved_psm_fdr * 100.0) << "%";
        }
        else
        {
          os << " (0 PSMs retained)";
        }
      }
      else
      {
        os << "  | PSM FDR filtering: off";
      }
      os << "\n";
    }

    // Charge distribution
    if (!s.charge_histogram.empty())
    {
      os << "[ProSE] Charges      :";
      for (const auto& [z, c] : s.charge_histogram) { os << " " << z << ":" << c; }
      os << "\n";
    }
    // Missed cleavages
    if (!s.missed_cleavage_histogram.empty())
    {
      os << "[ProSE] Missed clv   :";
      for (const auto& [mc, c] : s.missed_cleavage_histogram) { os << " " << mc << ":" << c; }
      os << "\n";
    }
    // Score distribution
    if (s.score_stats_valid)
    {
      os << "[ProSE] HyperScore   : min " << std::fixed << std::setprecision(1) << s.hyperscore_min
         << "  median " << s.hyperscore_median << "  max " << s.hyperscore_max << "\n";
    }

    // -- Tolerance estimation --
    if (s.prec_tol_valid || s.frag_tol_valid)
    {
      os << "[ProSE] --------------------- Tolerance estimate ---------------------\n";
      if (s.prec_tol_valid)
      {
        os << "[ProSE] Precursor    : median " << std::fixed << std::setprecision(2) << s.prec_err_median
           << " ppm, MAD " << s.prec_err_mad << "  -> recommend " << static_cast<int>(s.prec_err_recommended) << " ppm\n";
      }
      if (s.frag_tol_valid)
      {
        os << "[ProSE] Fragment     : median " << std::fixed << std::setprecision(2) << s.frag_err_median << " ppm, MAD " << s.frag_err_mad;
        if (s.frag_err_after_curve) { os << " (after the curve; the errors as measured: calibration:fragment_error_offset/sigma)\n"; }
        else { os << "  -> recommend " << static_cast<int>(s.frag_err_recommended) << " ppm\n"; }
      }
    }

    // -- Timing --
    os << "[ProSE] ------------------------- Timing -----------------------------\n";
    os << "[ProSE] " << std::fixed << std::setprecision(1);
    if (s.seconds_calibration > 0.0) { os << "calib " << s.seconds_calibration << "s | "; }
    if (s.seconds_search > 0.0)      { os << "scoring " << s.seconds_search << "s | "; }
    if (s.seconds_fdr > 0.0)         { os << "fdr " << s.seconds_fdr << "s | "; }
    os << "index(shared) " << sh.seconds_index_build << "s\n";
    // the whole run (set by the caller): the search of the file (reading it, rescoring, FDR, ...) and the outputs
    if (sh.seconds_total > 0.0) { os << "[ProSE] " << renderTotalTime(sh) << "\n"; }

    // -- Modification discovery (open search only) --
    if (is_open_search)
    {
      renderModificationSummary(mod_analysis, os);
    }
  }
} // namespace OpenMS
