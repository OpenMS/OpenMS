// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include "ProSEImpl.h"

#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#endif
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  double ProSEAlgorithm::CandidatePoolStats_::zScore() const
  {
    if (count < 2) return 0.0;
    const double n = static_cast<double>(count);
    const double mean = sum / n;
    const double var = std::max(0.0, sumsq / n - mean * mean);
    if (var <= 0.0) return 0.0; // every candidate scored the same: the best one is not an outlier
    return (best - mean) / std::sqrt(var);
  }

  void ProSEAlgorithm::annotatePeptDeepFeatures_(const PeakMap& spectra,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    if (!peptdeep_enable_) { return; }
#ifdef WITH_ONNX
    startProgress(0, 1, "Adding PeptDeep rescoring features...");
    PeptDeepRescoring rescoring;
    Param p = rescoring.getParameters();
    // A bare name resolves against share/OpenMS, an absolute path is returned unchanged.
    // The hint is needed because the compiled-in data path is the *source* share/OpenMS, which has
    // no models/ -- they are downloaded into the build tree. On an installed OpenMS the hint is a
    // dead path and File::find falls through to the installed data path, where 'make install' puts
    // them. This is the same resolution PeptDeepLibraryPredictor uses, deliberately: not relative
    // to the executable, which would only hold for a tool sitting in bin/.
    const StringList model_dirs = {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"};
    p.setValue("ms2_model", File::find(peptdeep_ms2_model_, model_dirs));
    p.setValue("rt_model", File::find(peptdeep_rt_model_, model_dirs));
    p.setValue("instrument", peptdeep_instrument_);
    p.setValue("nce", peptdeep_nce_);
    p.setValue("rt_model_type", peptdeep_rt_model_type_);
    p.setValue("batch_size", peptdeep_batch_size_);
    // Inference otherwise runs at its own default while the search around it scales with
    // the thread count the user actually asked for.
#ifdef _OPENMP
    p.setValue("threads", std::max(1, omp_get_max_threads()));
#else
    p.setValue("threads", 1);
#endif
    rescoring.setParameters(p);
    rescoring.setLogType(getLogType());
    rescoring.annotate(spectra, protein_ids, peptide_ids);
    endProgress();
    // 'peptdeep:enable' is an explicit request: an input that ends up without a single
    // predicted PSM must not silently continue with search-engine features only. An input
    // without PSMs has nothing to predict: its (empty) result is reported as it is.
    if (rescoring.getNumPredictedPSMs() == 0 && rescoring.getNumPSMs() > 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "[ProSE] 'peptdeep:enable' is set, but PeptDeep prediction added features to 0 / "
        + std::to_string(rescoring.getNumPSMs()) + " PSMs. Check that 'annotate:PSM' includes fragment annotations "
        "('ALL' or 'fragment_annotation').");
    }
#else
    (void)spectra; (void)protein_ids; (void)peptide_ids;
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "[ProSE] '" + paramName_("peptdeep:enable") + "' is set, but this OpenMS was built without ONNX support "
      "(WITH_ONNX=OFF), so the PeptDeep prediction-based rescoring features cannot be computed. "
      "Use an OpenMS build with -DWITH_ONNX=ON or set '" + paramName_("peptdeep:enable") + "' to 'false'.");
#endif
  }

  namespace
  {
    /// Scratch buffers of alignAbsoluteTolerance_(), reused across the hits and spectra of a thread (rewritten by every call).
    struct AlignmentScratch_
    {
      std::vector<double> prev;         ///< scores of the computed cells of row i - 1
      std::vector<double> cur;          ///< scores of the computed cells of row i
      std::vector<unsigned char> dir;   ///< traceback direction of every computed cell, rows concatenated
      std::vector<Size> row_begin;      ///< first computed column of row i
      std::vector<Size> row_offset;     ///< position of row i in dir (entry i + 1 closes row i)
    };

    /**
      @brief Absolute-tolerance (Da) branch of SpectrumAlignment::getSpectrumAlignment on flat storage.

      Returns exactly the alignment of the banded dynamic programme in SpectrumAlignment.h, which keeps
      its cells in std::map<Size, std::map<Size, ...>> (two tree nodes and several lookups per cell).
      Identity: the reference computes, per row i, one contiguous run of columns (from the left border
      it carries along to the column where it leaves the band), in the same order and with the same
      expressions as below. Every other cell it reads is either a border cell (i * tolerance or
      j * tolerance, (0,0) = 0) or absent, in which case it substitutes (i + j) * tolerance: the same
      value, so one formula serves both. A traceback step onto an absent cell reads a default (0,0)
      entry there, which ends the walk.
    */
    void alignAbsoluteTolerance_(std::vector<std::pair<Size, Size>>& alignment,
                                 const MSSpectrum& s1, const MSSpectrum& s2,
                                 const double tolerance, AlignmentScratch_& scratch)
    {
      if (!s1.isSorted() || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
      }
      alignment.clear();

      const Size n1 = s1.size();
      const Size n2 = s2.size();
      // cell (i, j) outside the computed runs: border value or the reference's substitute for an absent cell
      const auto outside = [tolerance](const Size i_plus_j) { return i_plus_j == 0 ? 0.0 : i_plus_j * tolerance; };
      enum : unsigned char { DIAGONAL = 1, UP = 2, LEFT = 3 }; // predecessor (i-1,j-1), (i,j-1), (i-1,j)

      std::vector<double>& prev = scratch.prev;
      std::vector<double>& cur = scratch.cur;
      std::vector<unsigned char>& dir = scratch.dir;
      prev.clear();
      dir.clear();
      scratch.row_begin.assign(n1 + 2, 0);
      scratch.row_offset.assign(n1 + 2, 0);

      Size left_ptr(1);
      Size last_i(0), last_j(0);
      Size prev_begin(1), prev_end(0); // computed columns of row i - 1 (row 0: none)
      for (Size i = 1; i <= n1; ++i)
      {
        const double pos1(s1[i - 1].getMZ());
        const Size cur_begin = left_ptr;
        cur.clear();
        scratch.row_begin[i] = cur_begin;
        scratch.row_offset[i] = dir.size();

        for (Size j = cur_begin; j <= n2; ++j)
        {
          bool off_band(false);
          const double pos2(s2[j - 1].getMZ());
          const double diff_align = fabs(pos1 - pos2);

          // running off the right border of the band?
          if (pos2 > pos1 && diff_align > tolerance)
          {
            if (i < n1 && j < n2 && s1[i].getMZ() < pos2)
            {
              off_band = true;
            }
          }

          // can we tighten the left border of the band?
          if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
          {
            ++left_ptr;
          }

          double score_align = diff_align;
          if (j - 1 >= prev_begin && j - 1 <= prev_end)
          {
            score_align += prev[j - 1 - prev_begin];
          }
          else
          {
            score_align += outside(i - 1 + j - 1);
          }

          double score_up = tolerance;
          if (j > cur_begin)
          {
            score_up += cur.back();
          }
          else
          {
            score_up += outside(i + j - 1);
          }

          double score_left = tolerance;
          if (j >= prev_begin && j <= prev_end)
          {
            score_left += prev[j - prev_begin];
          }
          else
          {
            score_left += outside(i - 1 + j);
          }

          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
            cur.push_back(score_align);
            dir.push_back(DIAGONAL);
            last_i = i;
            last_j = j;
          }
          else if (score_up <= score_left)
          {
            cur.push_back(score_up);
            dir.push_back(UP);
          }
          else
          {
            cur.push_back(score_left);
            dir.push_back(LEFT);
          }

          if (off_band)
          {
            break;
          }
        }
        prev_begin = cur_begin;
        prev_end = cur_begin + cur.size() - 1;
        prev.swap(cur);
      }
      scratch.row_offset[n1 + 1] = dir.size();

      // do traceback
      Size i = last_i;
      Size j = last_j;
      while (i >= 1 && j >= 1)
      {
        const Size begin = scratch.row_begin[i];
        const Size width = scratch.row_offset[i + 1] - scratch.row_offset[i];
        if (j < begin || j >= begin + width)
        {
          // absent cell: the reference reads (0,0), which counts as a diagonal step only at (1,1)
          if (i == 1 && j == 1) alignment.emplace_back(0, 0);
          break;
        }
        const unsigned char d = dir[scratch.row_offset[i] + (j - begin)];
        if (d == DIAGONAL)
        {
          alignment.emplace_back(i - 1, j - 1);
          --i;
          --j;
        }
        else if (d == UP)
        {
          --j;
        }
        else
        {
          --i;
        }
      }
      std::reverse(alignment.begin(), alignment.end());
    }
  }

  void ProSEAlgorithm::postProcessHits_(const PeakMap& exp,
        std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
        const std::vector<CandidatePoolStats_>& pool_stats,
        std::vector<ProteinIdentification>& protein_ids,
        PeptideIdentificationList& peptide_ids,
        Size top_hits,
  //      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
  //      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
  //      Size max_variable_mods_per_peptide, TODO: what about this parameter?
        const StringList& modifications_fixed,
        const StringList& modifications_variable,
        Int peptide_missed_cleavages,
        double precursor_mass_tolerance,
        double fragment_mass_tolerance,
        const std::string& precursor_mass_tolerance_unit_ppm,
        const std::string& fragment_mass_tolerance_unit_ppm,
        const Int precursor_min_charge,
        const Int precursor_max_charge,
        const std::string& enzyme,
        const std::string& database_name,
        const PeakMap* evidence_spectra,
        const std::vector<Size>* evidence_parent,
        bool defer_peak_annotations,
        const MzRecalibration_* fragment_recalibration) const
  {
    // Candidate-pool features (delta score, z-score, candidate count) are derived
    // from @p pool_stats rather than from @p annotated_hits: scoreSpectraAgainstIndex_
    // already pruned each spectrum to max(top_hits, 2) candidates, so the hits left
    // here are the report set, not the search space. pool_stats summarises every
    // candidate that was scored, zero-scoring ones included, and is accumulated
    // across chunks by the chunked search paths.
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }

    // One value per spectrum, kept in side vectors so scoring does not have to carry
    // them on every AnnotatedHit_ candidate.
    std::vector<float> delta_scores(annotated_hits.size(), 0.0f);
    std::vector<float> hyperscore_zscores(annotated_hits.size(), 0.0f);
    std::vector<float> ln_num_candidates(annotated_hits.size(), 0.0f);
#pragma omp parallel for default(none) shared(annotated_hits, top_hits, pool_stats, delta_scores, hyperscore_zscores, ln_num_candidates)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      const CandidatePoolStats_& stats = pool_stats[scan_index];
      delta_scores[scan_index] = static_cast<float>(stats.deltaScore());
      hyperscore_zscores[scan_index] = static_cast<float>(stats.zScore());
      // log1p rather than log: a spectrum can end up with zero candidates, and
      // ln(0) is not representable. The +1 offset is order-preserving, so this
      // stays a monotone proxy for the size of the scored search space.
      ln_num_candidates[scan_index] = static_cast<float>(std::log1p(static_cast<double>(stats.count)));

      // sort and keep n best elements according to score (unchanged)
      auto& hits = annotated_hits[scan_index];
      Size topn = top_hits > hits.size() ? hits.size() : top_hits;
      std::partial_sort(hits.begin(), hits.begin() + topn, hits.end(), AnnotatedHit_::hasBetterScore);
      hits.resize(topn);
      hits.shrink_to_fit();
    }

    bool annotation_precursor_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_fragment_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_prefix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_suffix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_num_matched_peaks = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::NUM_MATCHED_PEAKS) != annotate_psm_.end();
    bool annotation_matched_prefix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS) != annotate_psm_.end();
    bool annotation_matched_suffix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS) != annotate_psm_.end();
    bool annotation_longest_ion_run = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE) != annotate_psm_.end();
    bool annotation_matched_ion_current = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT) != annotate_psm_.end();
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end();
    bool annotation_hyperscore_zscore = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::HYPERSCORE_ZSCORE) != annotate_psm_.end();
    bool annotation_ln_num_candidates = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LN_NUM_CANDIDATES) != annotate_psm_.end();
    bool annotation_matched_ion_current_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT_FRACTION) != annotate_psm_.end();
    bool annotation_complementary_ions_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::COMPLEMENTARY_IONS_FRACTION) != annotate_psm_.end();
    const bool annotation_local_evidence = param_.getValue("annotate:local_fragment_evidence").toBool();
    const bool evidence_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    const Size n_evidence = evidence_spectra != nullptr ? evidence_spectra->size() : 0;
    if (annotation_local_evidence && (evidence_spectra == nullptr
        || n_evidence + (evidence_parent != nullptr ? evidence_parent->size() : 0) != exp.size()))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Local fragment evidence requires the aligned peak lists retained during preprocessing.");
    }

    // "ALL" adds all annotations
    if (std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end())
    {
      annotation_precursor_error_ppm = true;
      annotation_fragment_error_ppm = true;
      annotation_prefix_fraction = true;
      annotation_suffix_fraction = true;
      annotation_num_matched_peaks = true;
      annotation_matched_prefix_ions = true;
      annotation_matched_suffix_ions = true;
      annotation_longest_ion_run = true;
      annotation_matched_ion_current = true;
      annotation_fragment_annotations = true;
      annotation_hyperscore_zscore = true;
      annotation_ln_num_candidates = true;
      annotation_matched_ion_current_fraction = true;
      annotation_complementary_ions_fraction = true;
    }

    // Alignment is needed for fragment error, fragment annotations, longest ion run, MIC,
    // normalized MIC, and complementary ion pairs
    // (deferred fragment annotations are aligned when they are made)
    const bool need_alignment = annotation_fragment_error_ppm || (annotation_fragment_annotations && !defer_peak_annotations)
      || annotation_longest_ion_run || annotation_matched_ion_current || annotation_matched_ion_current_fraction
      || annotation_complementary_ions_fraction;
    // no ion name is written here: the theoretical spectra may come without them
    const bool without_ion_names = !annotation_fragment_annotations || defer_peak_annotations;

    // Both configurations depend only on the fragment tolerance, so they are
    // shared read-only by every thread of the loop below: getSpectrum() and
    // getSpectrumAlignment() are const and hold no mutable state (the scoring
    // loop in scoreSpectraAgainstIndex_ shares its generators the same way).
    // The alignment tolerance mirrors the search's fragment tolerance so the
    // reported FRAGMENT_ERROR_MEDIAN_PPM is not polluted by far-off spurious
    // matches — SpectrumAlignment's default is 0.3 Da absolute, which is ~30×
    // looser than a typical 20 ppm search.
    // The generators mirror the ion series of the scoring, so annotation features (fragment
    // error, peak annotations, longest ion run, MIC) are computed on the same theoretical
    // spectrum the candidate was scored against.
    const SpectrumGenerators_ generators = spectrumGenerators_();
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }
    // Da tolerance: same alignment from a flat-storage copy of SpectrumAlignment's banded DP
    // (alignAbsoluteTolerance_ above); ppm tolerance keeps SpectrumAlignment's cheap matching.
    const bool sa_absolute = !sa.getParameters().getValue("is_relative_tolerance").toBool();
    const double sa_tolerance = (double)sa.getParameters().getValue("tolerance");

    // Meta value keys of the loop below, resolved to registry indices once per call: every
    // string-keyed setMetaValue()/getMetaValue() takes the registry's process-wide lock (also for
    // names that are already registered), which serialised the annotation threads.
    // Identity: MetaInfo keeps (and idXML writes) the values of a hit ordered by registry index, so
    // new names must get their indices in the order in which the loop registered them: the first
    // spectrum with hits registers "scan_index", then IM if it carries a drift time, then its first
    // hit registers the per-hit names in the order they are set below. All hits set the same names,
    // so no thread could register a later name before an earlier one and this order did not depend
    // on the thread count. registerName() only looks up names that are already known, and without
    // any hit nothing was registered, so nothing is registered here either. Names that stay
    // string-keyed (spectrum reference, IM, rank; at most one lookup per spectrum or hit) are
    // still registered inside the loop, after the ones below, as before (checked: the registry
    // contents after a search are the same as before this change, at 1 and at 16 threads).
    const int isotope_error_sign = isotope_error_observed_minus_theoretical_ ? -1 : 1;
    const bool open_search_mode = isOpenSearchMode_();
    // calibration:recalibrate kept the observed precursor m/z under this name (registered by applyRecalibration_(); UInt(-1):
    // never applied, and the registry stays as it is)
    const UInt mv_observed_mz = MetaInfoInterface::metaRegistry().getIndex(OBSERVED_PRECURSOR_MZ);
    UInt mv_scan_index{}, mv_fragment_error{}, mv_precursor_error{}, mv_prefix_fraction{}, mv_suffix_fraction{},
         mv_num_matched_peaks{}, mv_matched_prefix_ions{}, mv_matched_suffix_ions{}, mv_delta_score{},
         mv_hyperscore_zscore{}, mv_ln_num_candidates{}, mv_matched_ion_current{}, mv_matched_ion_current_fraction{},
         mv_longest_ion_run{}, mv_complementary_ions_fraction{}, mv_chance_match_surprise{},
         mv_mass_competition_evidence{}, mv_isotope_error{}, mv_delta_mass{}, mv_rank{};
    {
      const auto first_with_hits = std::find_if(annotated_hits.begin(), annotated_hits.end(),
                                                [](const std::vector<AnnotatedHit_>& hits) { return !hits.empty(); });
      if (first_with_hits != annotated_hits.end())
      {
        MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
        mv_scan_index = registry.registerName("scan_index");
        const MSSpectrum& first_spectrum = exp[first_with_hits - annotated_hits.begin()];
        if (IMTypes::determineIMFormat(first_spectrum) == IMFormat::IM_SPECTRUM)
        {
          registry.registerName(first_spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE
                                  ? Constants::UserParam::FAIMS_CV : Constants::UserParam::IM);
        }
        if (annotation_fragment_error_ppm) mv_fragment_error = registry.registerName(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
        if (annotation_precursor_error_ppm) mv_precursor_error = registry.registerName(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
        if (annotation_prefix_fraction) mv_prefix_fraction = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
        if (annotation_suffix_fraction) mv_suffix_fraction = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
        if (annotation_num_matched_peaks) mv_num_matched_peaks = registry.registerName(Constants::UserParam::NUM_MATCHED_PEAKS);
        if (annotation_matched_prefix_ions) mv_matched_prefix_ions = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS);
        if (annotation_matched_suffix_ions) mv_matched_suffix_ions = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS);
        mv_delta_score = registry.registerName(Constants::UserParam::DELTA_SCORE);
        if (annotation_hyperscore_zscore) mv_hyperscore_zscore = registry.registerName(Constants::UserParam::HYPERSCORE_ZSCORE);
        if (annotation_ln_num_candidates) mv_ln_num_candidates = registry.registerName(Constants::UserParam::LN_NUM_CANDIDATES);
        if (annotation_matched_ion_current) mv_matched_ion_current = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT);
        if (annotation_matched_ion_current_fraction) mv_matched_ion_current_fraction = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
        if (annotation_longest_ion_run) mv_longest_ion_run = registry.registerName(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
        if (annotation_complementary_ions_fraction) mv_complementary_ions_fraction = registry.registerName(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
        if (annotation_local_evidence)
        {
          mv_chance_match_surprise = registry.registerName(Constants::UserParam::CHANCE_MATCH_SURPRISE);
          mv_mass_competition_evidence = registry.registerName(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
        }
        mv_isotope_error = registry.registerName(Constants::UserParam::ISOTOPE_ERROR);
        if (open_search_mode) mv_delta_mass = registry.registerName("DeltaMass");
        // The rank of the hits below the top one (PeptideHit::setRank() by name: one lookup under the registry's
        // lock per hit), registered where the loop below registered it (by its first second-best hit), after the names
        // above: no other name of a hit is registered in between.
        if (std::any_of(annotated_hits.begin(), annotated_hits.end(), [](const std::vector<AnnotatedHit_>& hits) { return hits.size() > 1; }))
        {
          mv_rank = registry.registerName(Constants::UserParam::RANK);
        }
      }
    }

    // The meta values a hit gets below, and later in the search (protein mapping: 2, ion priors: 2): reserved at once.
    const Size hit_meta_values = Size(annotation_fragment_error_ppm) + annotation_precursor_error_ppm + annotation_prefix_fraction
      + annotation_suffix_fraction + annotation_num_matched_peaks + annotation_matched_prefix_ions + annotation_matched_suffix_ions
      + 1 /* delta score */ + annotation_hyperscore_zscore + annotation_ln_num_candidates + annotation_matched_ion_current
      + annotation_matched_ion_current_fraction + annotation_longest_ion_run + annotation_complementary_ions_fraction
      + 2 * Size(annotation_local_evidence) + 1 /* isotope error */ + open_search_mode + 1 /* rank */ + 4;

    // With an empty list to start from, every spectrum with hits gets a slot in scan order: the identifications
    // end up in scan order without a lock and without the sort below (which orders them by scan_index, i.e. the
    // same order). Otherwise they are appended under a lock and sorted, as before.
    const bool slots = peptide_ids.empty();
    std::vector<Size> slot_of_scan;
    if (slots)
    {
      slot_of_scan.resize(annotated_hits.size());
      Size used_slots = 0;
      for (Size scan_index = 0; scan_index < annotated_hits.size(); ++scan_index)
      {
        slot_of_scan[scan_index] = used_slots;
        if (!annotated_hits[scan_index].empty()) ++used_slots;
      }
      peptide_ids.resize(used_slots);
    }

    // The body of the parallel region below. A lambda, so that the region captures one variable instead of
    // dozens: LLVM's OpenMP runtime passes at most 32 to a region on WebAssembly. Its loop is an orphaned
    // worksharing construct, bound to the region that calls the lambda.
    const auto annotate_scans = [&]()
    {
    // scratch of the theoretical spectra, alignments and local evidence, reused by all hits of all spectra of a thread
    // (as before by all hits of a spectrum: each use overwrites what it reads)
    AlignmentScratch_ alignment_scratch;
    std::vector<double> evidence_alternatives;
    // the theoretical spectrum without ion names, see below
    MSSpectrum coded_spec;
    std::vector<uint8_t> coded_series;
    std::vector<UInt> coded_ordinals;
    std::vector<Int> coded_charges;
#pragma omp for schedule(dynamic, 16)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spec = exp[scan_index];
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setSpectrumReference( spec.getNativeID());
        pi.setMetaValue(mv_scan_index, static_cast<unsigned int>(scan_index));
        pi.setScoreType(SEARCH_ENGINE_SCORE); // also for scoring:method mass_accuracy, a log-space HyperScore (recorded in the search parameters)
        pi.setHigherScoreBetter(true);
        const Precursor& precursor = spec.getPrecursors()[0];
        double mz = precursor.getMZ(); // the m/z the search used (corrected by calibration:recalibrate, if on)
        pi.setRT(spec.getRT());
        // the PSM reports the precursor m/z as observed
        pi.setMZ(mv_observed_mz != UInt(-1) && precursor.metaValueExists(mv_observed_mz) ? static_cast<double>(precursor.getMetaValue(mv_observed_mz)) : mz);

        // Annotate the spectrum's single drift time: a FAIMS compensation voltage as FAIMS_CV (which the
        // FAIMS-aware tools read), an ion mobility (e.g. DDA-PASEF, 1/K0) as IM
        if (IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM)
        {
          pi.setMetaValue(spec.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE
                            ? Constants::UserParam::FAIMS_CV : Constants::UserParam::IM,
                          spec.getDriftTime());
        }

        Size charge = spec.getPrecursors()[0].getCharge();

        // Spectrum-level quantity, identical for every candidate of this spectrum, so it is
        // computed once here rather than per hit.
        const double spectrum_tic =
          annotation_matched_ion_current_fraction ? spec.calculateTIC() : 0.0;
        const MSSpectrum& evidence_spec = !annotation_local_evidence ? spec
          : (*evidence_spectra)[static_cast<Size>(scan_index) < n_evidence ? static_cast<Size>(scan_index)
                                                                           : (*evidence_parent)[scan_index - n_evidence]];
        const std::vector<double> local_densities = annotation_local_evidence
          ? localPeakDensities_(evidence_spec) : std::vector<double>{};

        // create full peptide hit structure from annotated hits (their sequences are moved: the hits are freed after
        // this function)
        vector<PeptideHit> phs;
        phs.reserve(annotated_hits[scan_index].size());
        for (auto& ah : annotated_hits[scan_index])
        {
          PeptideHit ph;
          ph.reserveMetaValues(hit_meta_values); // (one block instead of growing by doubling)
          // Prefer spectrum charge; if absent (0), fall back to the charge actually used by FI for this candidate
          const Size used_charge = (charge > 0) ? charge : static_cast<Size>(ah.applied_charge);
          ph.setCharge(used_charge);
          ph.setScore(ah.score);
          ph.setSequence(std::move(ah.sequence));
          const AASequence& sequence = ph.getSequence();

          // Generate theoretical spectrum + alignment for annotations that need it. Unless the ion names are written
          // (fragment annotations here), the peaks of getSpectrum() come without them, with the series, ordinal and
          // charge of each ion as numbers (getIonSeriesSpectrum(): the same peaks in the same order; about a fifth of
          // the cost), if the generator supports that.
          std::vector<std::pair<Size, Size>> alignment;
          MSSpectrum theoretical_spec;
          bool coded = false;
          // Annotate the charges actually scored when higher charges are scored.
          const int max_frag_z = scoring_multiple_charges_ ? scoringMaxCharge_(static_cast<int>(used_charge))
                                                           : ((charge >= 2) ? std::min<int>(charge - 1, 2) : 1);
          if (need_alignment)
          {
            coded = without_ion_names
                    && tsg.getIonSeriesSpectrum(coded_spec, coded_series, coded_ordinals, coded_charges, sequence, 1, max_frag_z);
            if (!coded) tsg.getSpectrum(theoretical_spec, sequence, 1, max_frag_z);
            if (sa_absolute)
            {
              alignAbsoluteTolerance_(alignment, coded ? coded_spec : theoretical_spec, spec, sa_tolerance, alignment_scratch);
            }
            else
            {
              sa.getSpectrumAlignment(alignment, coded ? coded_spec : theoretical_spec, spec);
            }
          }
          const MSSpectrum& theo_spec = coded ? coded_spec : theoretical_spec;

          if (annotation_fragment_error_ppm)
          {
            std::vector<double> err;
            for (const auto& match : alignment)
            {
              double fragment_error = fabs(Math::getPPM(spec[match.second].getMZ(), theo_spec[match.first].getMZ()));
              err.push_back(fragment_error);
            }
            double median_ppm_error(0);
            if (!err.empty()) { median_ppm_error = Math::median(err.begin(), err.end(), false); }
            ph.setMetaValue(mv_fragment_error, median_ppm_error);
          }

          if (annotation_precursor_error_ppm)
          {
            // Subtract out the isotope offset FI matched at — FragmentIndex searches
            // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs
            // + isotope_error * C13C12, and the observed-to-monoiso correction in m/z is
            //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
            // (ah.isotope_error is this search offset, whatever sign the PSM reports).
            // Without this, a ±1 Da FI match reports ~1000 ppm / charge for the Percolator
            // feature, corrupting target/decoy discrimination.
            const double corrected_mz = monoisotopicPrecursorMZ(mz, ah.isotope_error, static_cast<double>(used_charge));
            double theo_mz = sequence.getMZ(used_charge);
            double ppm_difference = Math::getPPM(corrected_mz, theo_mz);
            ph.setMetaValue(mv_precursor_error, ppm_difference);
          }

          if (annotation_prefix_fraction)
          {
            ph.setMetaValue(mv_prefix_fraction, ah.prefix_fraction);
          }

          if (annotation_suffix_fraction)
          {
            ph.setMetaValue(mv_suffix_fraction, ah.suffix_fraction);
          }

          // Matched ion counts (from scoring, no alignment needed)
          if (annotation_num_matched_peaks)
          {
            ph.setMetaValue(mv_num_matched_peaks, static_cast<int>(ah.matched_prefix_ions + ah.matched_suffix_ions));
          }
          if (annotation_matched_prefix_ions)
          {
            ph.setMetaValue(mv_matched_prefix_ions, static_cast<int>(ah.matched_prefix_ions));
          }
          if (annotation_matched_suffix_ions)
          {
            ph.setMetaValue(mv_matched_suffix_ions, static_cast<int>(ah.matched_suffix_ions));
          }

          ph.setMetaValue(mv_delta_score, delta_scores[scan_index]);

          if (annotation_hyperscore_zscore)
          {
            ph.setMetaValue(mv_hyperscore_zscore, hyperscore_zscores[scan_index]);
          }
          if (annotation_ln_num_candidates)
          {
            ph.setMetaValue(mv_ln_num_candidates, ln_num_candidates[scan_index]);
          }

          // Fragment annotations, longest ion run, MIC, normalized MIC, and complementary
          // ion pairs all iterate the alignment + ion names
          if ((annotation_fragment_annotations && !defer_peak_annotations) || annotation_longest_ion_run || annotation_matched_ion_current
            || annotation_matched_ion_current_fraction || annotation_complementary_ions_fraction)
          {
            // (the names and charges of getSpectrum(); none for a coded spectrum)
            const auto* ion_names = theoretical_spec.getStringDataArrays().empty() ? nullptr : &theoretical_spec.getStringDataArrays()[0];
            const auto* ion_charges = theoretical_spec.getIntegerDataArrays().empty() ? nullptr : &theoretical_spec.getIntegerDataArrays()[0];

            // Build PeakAnnotation vector + collect ion ordinals for longest run.
            // Prefix = a/b/c (N-terminal), suffix = x/y/z (C-terminal). Ordinals
            // for different ion types at the same cleavage position (e.g. a3 + b3)
            // are merged via std::unique below — each ordinal is a backbone
            // position, not an ion-type-specific identifier.
            std::vector<PeptideHit::PeakAnnotation> peak_annotations;
            std::vector<int> prefix_ordinals, suffix_ordinals;
            double matched_ion_current = 0.0;
            const bool need_mic = annotation_matched_ion_current || annotation_matched_ion_current_fraction;
            const bool need_ordinals = annotation_longest_ion_run || annotation_complementary_ions_fraction;
            // Dedup guard for MIC: in ppm-alignment mode a single experimental
            // peak can match multiple theoretical peaks (e.g. b-ion and near
            // isotope), so we must sum each exp_idx at most once. Sized only
            // when MIC is actually requested.
            std::vector<char> counted_exp_peaks(need_mic ? spec.size() : 0, 0);
            if (annotation_fragment_annotations && !defer_peak_annotations) { peak_annotations.reserve(alignment.size()); }

            for (const auto& [theo_idx, exp_idx] : alignment)
            {
              if (annotation_fragment_annotations && !defer_peak_annotations)
              {
                PeptideHit::PeakAnnotation pa;
                // the m/z as observed (calibration:recalibrate corrected the peak in place), as the PSM's precursor m/z
                pa.mz = fragment_recalibration != nullptr ? fragment_recalibration->uncorrect(spec[exp_idx].getMZ(), spec.getRT())
                                                          : spec[exp_idx].getMZ();
                pa.intensity = spec[exp_idx].getIntensity();
                pa.annotation = (*ion_names)[theo_idx];
                pa.charge = (*ion_charges)[theo_idx];
                peak_annotations.push_back(pa);
              }

              if (need_mic && !counted_exp_peaks[exp_idx])
              {
                matched_ion_current += spec[exp_idx].getIntensity();
                counted_exp_peaks[exp_idx] = 1;
              }

              if (need_ordinals && coded)
              {
                // the number of the ion name ("b5+", "z.4+": 5, 4), never 0
                (coded_series[theo_idx] == 0 ? prefix_ordinals : suffix_ordinals).push_back(static_cast<int>(coded_ordinals[theo_idx]));
              }
              else if (need_ordinals && (*ion_names)[theo_idx].size() >= 2)
              {
                const std::string& name = (*ion_names)[theo_idx];
                const char c = name[0];
                const bool is_prefix = (c == 'a' || c == 'b' || c == 'c');
                const bool is_suffix = (c == 'x' || c == 'y' || c == 'z');
                if (is_prefix || is_suffix)
                {
                  // Extract ordinal: "b5", "y3-H2O1+", "c12++", "z.4+" (z+1) -> 5, 3, 12, 4
                  Size pos = 1;
                  while (pos < name.size() && (name[pos] == '.' || name[pos] == '\'')) ++pos; // z. (z+1), z' (z+2)
                  const Size ordinal_begin = pos;
                  while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
                  if (pos > ordinal_begin)
                  {
                    int ordinal = StringUtils::toInt32(StringUtils::substr(name, ordinal_begin, pos - ordinal_begin));
                    (is_prefix ? prefix_ordinals : suffix_ordinals).push_back(ordinal);
                  }
                }
              }
            }

            if (annotation_fragment_annotations && !defer_peak_annotations)
            {
              ph.setPeakAnnotations(std::move(peak_annotations));
            }

            if (annotation_matched_ion_current)
            {
              ph.setMetaValue(mv_matched_ion_current, matched_ion_current);
            }

            if (annotation_matched_ion_current_fraction)
            {
              ph.setMetaValue(mv_matched_ion_current_fraction,
                              spectrum_tic > 0 ? matched_ion_current / spectrum_tic : 0.0);
            }

            if (need_ordinals)
            {
              // Compute longest consecutive run across prefix and suffix series.
              // Sorts + deduplicates each ordinal vector in place (a backbone position
              // matched by multiple ion types, e.g. a3 and b3, counts once).
              auto longestRun = [](std::vector<int>& v) -> int {
                if (v.empty()) return 0;
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                int best = 1, run = 1;
                for (Size i = 1; i < v.size(); ++i)
                {
                  if (v[i] == v[i - 1] + 1) { ++run; if (run > best) best = run; }
                  else run = 1;
                }
                return best;
              };
              int longest_prefix = longestRun(prefix_ordinals);
              int longest_suffix = longestRun(suffix_ordinals);

              if (annotation_longest_ion_run)
              {
                ph.setMetaValue(mv_longest_ion_run, std::max(longest_prefix, longest_suffix));
              }

              if (annotation_complementary_ions_fraction)
              {
                // A prefix ion at backbone position i (a_i/b_i/c_i) is complementary to a
                // suffix ion at position (peptide_length - i) (x/y/z at the same cleavage
                // site). Andromeda-style structural signal, distinct from the independently
                // computed prefix/suffix fractions above.
                const int pep_len = static_cast<int>(sequence.size());
                double complementary_fraction = 0.0;
                if (pep_len > 1)
                {
                  Size n_complementary = 0;
                  for (int p : prefix_ordinals)
                  {
                    if (std::binary_search(suffix_ordinals.begin(), suffix_ordinals.end(), pep_len - p)) { ++n_complementary; }
                  }
                  complementary_fraction = static_cast<double>(n_complementary) / static_cast<double>(pep_len - 1);
                }
                ph.setMetaValue(mv_complementary_ions_fraction, complementary_fraction);
              }
            }
          }


          if (annotation_local_evidence)
          {
            // Use the fragment charges that were scored. The annotation spectrum (charges 1..max_frag_z)
            // contains them whenever it was generated, so it is reused instead of generating a second one.
            const int evidence_z = scoringMaxCharge_(static_cast<int>(used_charge));
            LocalFragmentEvidence_ evidence;
            if (need_alignment && evidence_z <= max_frag_z)
            {
              evidence = coded ? localFragmentEvidence_(evidence_spec, coded_spec, coded_charges, evidence_z, local_densities,
                                                        fragment_mass_tolerance, evidence_ppm, evidence_alternatives)
                               : localFragmentEvidence_(evidence_spec, theoretical_spec, evidence_z, local_densities,
                                                        fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            else if (tsg.getIonSeriesSpectrum(coded_spec, coded_series, coded_ordinals, coded_charges, sequence, 1, evidence_z))
            {
              // (the peaks and charges of getSpectrum(), without the names)
              evidence = localFragmentEvidence_(evidence_spec, coded_spec, coded_charges, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            else
            {
              MSSpectrum evidence_theory;
              tsg.getSpectrum(evidence_theory, sequence, 1, evidence_z);
              evidence = localFragmentEvidence_(evidence_spec, evidence_theory, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            ph.setMetaValue(mv_chance_match_surprise, evidence.chance_match_surprise);
            ph.setMetaValue(mv_mass_competition_evidence, evidence.mass_competition_evidence);
          }

          // Add isotope error metavalue (always; exposed as Percolator feature). ah.isotope_error is the offset
          // FragmentIndex added to the observed mass (theoretical minus observed); report:isotope_error_convention
          // selects the reported sign.
          ph.setMetaValue(mv_isotope_error, isotope_error_sign * ah.isotope_error);

          // Add delta mass metavalue for open search
          if (open_search_mode)
          {
            ph.setMetaValue(mv_delta_mass, ah.delta_mass);
          }

          // store PSM
          phs.push_back(std::move(ph));
        }
        pi.setHits(std::move(phs));
        // Ensure hits are sorted by score (best first), then assign ranks explicitly (0 = top hit)
        pi.sort();
        {
          std::vector<PeptideHit>& hits = pi.getHits();
          for (Size r = 0; r < hits.size(); ++r)
          {
            // as setRank(r): the rank of the top hit (0) is the default, which a new hit does not carry
            if (r != 0) hits[r].setMetaValue(mv_rank, static_cast<UInt>(r));
          }
        }

        // Debug: log spectrum-level top hit details before storing PeptideIdentification.
        // DEBUG (not INFO) because multi-file mode would emit one line per scan per input.
        if (!pi.getHits().empty())
        {
          const PeptideHit& top_hit = pi.getHits().front();
          OPENMS_LOG_DEBUG << "[ProSE] scan_index=" << scan_index
                           << " top_ln(hyperscore)=" << top_hit.getScore()
                           << " top_charge=" << top_hit.getCharge()
                           << " top_isotope_error=" << (int)top_hit.getMetaValue(mv_isotope_error)
                           << std::endl;
        }
        if (slots)
        {
          peptide_ids[slot_of_scan[scan_index]] = std::move(pi);
        }
        else
        {
#pragma omp critical (peptide_ids_access)
          {
            //clang-tidy: seems to be a false-positive in combination with omp
            peptide_ids.push_back(std::move(pi));
          }
        }
      }
    }
    };

    // The work per spectrum varies with its number of hits: hand out small blocks
#pragma omp parallel
    {
      annotate_scans();
    } // omp parallel

#ifdef _OPENMP
    // we need to sort the peptide_ids by scan_index in order to have the same output in the idXML-file
    if (omp_get_max_threads() > 1 && !slots)
    {
      // one registry lookup for the whole sort instead of two (locked) lookups per comparison;
      // getMetaValue(name) is getMetaValue(getIndex(name)), so the comparisons are unchanged
      const UInt scan_index_key = MetaInfoInterface::metaRegistry().getIndex("scan_index");
      std::sort(peptide_ids.begin(), peptide_ids.end(), [scan_index_key](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return a.getMetaValue(scan_index_key) < b.getMetaValue(scan_index_key);
      });
    }
#endif

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("ProSE");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    DateTime now = DateTime::now();
    std::string identifier("ProSE_" + now.get());
    protein_ids[0].setIdentifier(identifier);
    for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges =StringUtils::toStr(precursor_min_charge) + ":" + StringUtils::toStr(precursor_max_charge);

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = modifications_fixed;
    search_parameters.variable_modifications = modifications_variable;
    search_parameters.missed_cleavages = peptide_missed_cleavages;
    search_parameters.fragment_mass_tolerance = fragment_mass_tolerance;
    search_parameters.precursor_mass_tolerance = precursor_mass_tolerance;
    search_parameters.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);

    // add additional percolator features or post-processing
    StringList feature_set{"score"};
    if (annotation_fragment_error_ppm) feature_set.push_back(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
    if (annotation_prefix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
    if (annotation_suffix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
    if (annotation_longest_ion_run) feature_set.push_back(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
    if (annotation_matched_prefix_ions) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS);
    if (annotation_matched_suffix_ions) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS);
    if (annotation_matched_ion_current) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT);
    if (annotation_matched_ion_current_fraction) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
    if (annotation_complementary_ions_fraction) feature_set.push_back(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
    if (annotation_hyperscore_zscore) feature_set.push_back(Constants::UserParam::HYPERSCORE_ZSCORE);
    if (annotation_ln_num_candidates) feature_set.push_back(Constants::UserParam::LN_NUM_CANDIDATES);
    if (annotation_local_evidence)
    {
      feature_set.push_back(Constants::UserParam::CHANCE_MATCH_SURPRISE);
      feature_set.push_back(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
    }
    feature_set.push_back(Constants::UserParam::DELTA_SCORE);
    feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
    if (chimericSearch_()) feature_set.push_back(CHIMERIC_RANK); // set by search() after this
    // note: precursor error is calculated by percolator itself
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    // Readers tell the sign of the PSMs' isotope_error by this record; without it (earlier ProSE versions, or
    // report:isotope_error_convention=theoretical_minus_observed) the sign is theoretical minus observed.
    if (isotope_error_observed_minus_theoretical_)
    {
      search_parameters.setMetaValue("isotope_error_convention", "observed_minus_theoretical");
    }
    // record whether open-search mode was used
    search_parameters.setMetaValue("open_search", isOpenSearchMode_() ? "true" : "false");

    search_parameters.setMetaValue("peptide:clip_nterm_methionine", param_.getValue("peptide:clip_nterm_methionine"));
    search_parameters.setMetaValue("peptide:deduplicate", param_.getValue("peptide:deduplicate"));
    // recorded only when not the default (the default output stays as before)
    if (decoy_keep_termini_) search_parameters.setMetaValue("decoy_method", "keep_termini");
    if (nterm_options_disabled_) // the effective values differ from the requested ones (buildDecoyAugmentedDB_)
    {
      search_parameters.setMetaValue("rescoring:terminal_features", "false");
      search_parameters.setMetaValue("rescoring:stratum_q_values", "false");
      search_parameters.setMetaValue("rescoring:nterm_options_disabled", "the database's decoys rarely start with M");
    }
    search_parameters.setMetaValue("fragment:query_spectrum", param_.getValue("fragment:query_spectrum"));
    search_parameters.setMetaValue("fragment:query_spectrum_resolved", query_raw_spectrum_ ? "raw" : "processed");
    search_parameters.setMetaValue("annotate:local_fragment_evidence", param_.getValue("annotate:local_fragment_evidence"));
    search_parameters.setMetaValue("scoring:fragment_charges", param_.getValue("scoring:fragment_charges"));
    search_parameters.setMetaValue("scoring:fragment_charges_resolved", scoring_multiple_charges_ ? "multiple" : "single");
    if (mass_accuracy_score_) // recorded when it changes the score, so HyperScore searches keep their parameter list
    {
      search_parameters.setMetaValue("scoring:method", param_.getValue("scoring:method"));
      search_parameters.setMetaValue("scoring:method_resolved", "mass_accuracy");
      search_parameters.setMetaValue("scoring:mass_error_sd", mass_error_sd_ppm_);
    }
    search_parameters.setMetaValue("peaks:window_type", peaks_window_type_);
    search_parameters.setMetaValue(
      "peaks:window_type_resolved",
      peaks_window_type_ == "auto"
        ? (Deisotoper::isToleranceSupported(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm") ? "jump_full" : "jump")
        : peaks_window_type_);
    // MS2 deisotoping rule (used only when the spectra are deisotoped, see fragment:deisotope)
    search_parameters.setMetaValue("fragment:deisotope_min_peaks", static_cast<int>(deisotoping_.min_peaks));
    search_parameters.setMetaValue("fragment:deisotope_charge_cap", deisotoping_.charge_cap_precursor ? "precursor" : "none");
    search_parameters.setMetaValue("fragment:deisotope_sum_intensity", deisotoping_.sum_intensity ? "true" : "false");

    search_parameters.enzyme_term_specificity = peptide_enzyme_specificity_;
    protein_ids[0].setSearchParameters(std::move(search_parameters));

    // Annotate the IM unit on the ProteinIdentification if all PeptideIdentifications have an IM value. The
    // unit is that of the spectra the values came from (FAIMS compensation voltages are written as FAIMS_CV,
    // so they never count here)
    if (!peptide_ids.empty())
    {
      bool all_have_im = std::all_of(
        peptide_ids.begin(), peptide_ids.end(),
        [](const PeptideIdentification& pid) { return pid.metaValueExists(Constants::UserParam::IM); });

      if (all_have_im)
      {
        const auto with_im = std::find_if(exp.begin(), exp.end(), [](const MSSpectrum& spec) {
          return IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM
                 && spec.getDriftTimeUnit() != DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
        });
        protein_ids[0].setMetaValue(Constants::UserParam::IM,
          (with_im != exp.end() ? *with_im : exp[0]).getDriftTimeUnitAsString());
      }
    }
  }

  // =====================================================================
  // Score all spectra against a single FragmentIndex, appending results to
  // annotated_hits. Used by both chunked and non-chunked search paths.
  // =====================================================================
  void ProSEAlgorithm::annotateFragmentPeaks_(const PeakMap& exp, PeptideIdentificationList& peptide_ids,
                                              double fragment_mass_tolerance,
                                              const std::string& fragment_mass_tolerance_unit_ppm,
                                              const MzRecalibration_* fragment_recalibration) const
  {
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end()
      || std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end();
    if (!annotation_fragment_annotations) return;
    // the generators and the alignment of postProcessHits_()
    const SpectrumGenerators_ generators = spectrumGenerators_();
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }
    const bool sa_absolute = !sa.getParameters().getValue("is_relative_tolerance").toBool();
    const double sa_tolerance = (double)sa.getParameters().getValue("tolerance");
    const UInt scan_index_key = MetaInfoInterface::metaRegistry().getIndex("scan_index");

    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      try
      {
        PeptideIdentification& pi = peptide_ids[i];
        if (pi.getHits().empty()) continue;
        const Size scan_index = static_cast<Size>(static_cast<long long>(pi.getMetaValue(scan_index_key))); // (set by postProcessHits_())
        const MSSpectrum& spec = exp[scan_index];
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);
        const Size charge = spec.getPrecursors()[0].getCharge();
        AlignmentScratch_ alignment_scratch;
        for (PeptideHit& ph : pi.getHits())
        {
          // as postProcessHits_() for this hit (its charge is the charge used there)
          const Size used_charge = static_cast<Size>(ph.getCharge());
          const int max_frag_z = scoring_multiple_charges_ ? scoringMaxCharge_(static_cast<int>(used_charge))
                                                           : ((charge >= 2) ? std::min<int>(charge - 1, 2) : 1);
          MSSpectrum theoretical_spec;
          std::vector<std::pair<Size, Size>> alignment;
          tsg.getSpectrum(theoretical_spec, ph.getSequence(), 1, max_frag_z);
          if (sa_absolute)
          {
            alignAbsoluteTolerance_(alignment, theoretical_spec, spec, sa_tolerance, alignment_scratch);
          }
          else
          {
            sa.getSpectrumAlignment(alignment, theoretical_spec, spec);
          }
          const auto& ion_names = theoretical_spec.getStringDataArrays()[0];
          const auto& ion_charges = theoretical_spec.getIntegerDataArrays()[0];
          std::vector<PeptideHit::PeakAnnotation> peak_annotations;
          peak_annotations.reserve(alignment.size());
          for (const auto& [theo_idx, exp_idx] : alignment)
          {
            PeptideHit::PeakAnnotation pa;
            // the m/z as observed, as in postProcessHits_()
            pa.mz = fragment_recalibration != nullptr ? fragment_recalibration->uncorrect(spec[exp_idx].getMZ(), spec.getRT())
                                                      : spec[exp_idx].getMZ();
            pa.intensity = spec[exp_idx].getIntensity();
            pa.annotation = ion_names[theo_idx];
            pa.charge = ion_charges[theo_idx];
            peak_annotations.push_back(pa);
          }
          ph.setPeakAnnotations(std::move(peak_annotations));
        }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
  }

  void ProSEAlgorithm::annotateIonPriors_(const PeakMap& spectra,
                                          const FragmentIonLikelihoodModel::PeakLists& evidence,
                                          std::vector<ProteinIdentification>& protein_ids,
                                          PeptideIdentificationList& peptide_ids) const
  {
    if (!self_trained_ion_priors_) return;
    StopWatch sw_train, sw_annotate;
    sw_train.start();
    const bool ppm = fragment_mass_tolerance_unit_ == "ppm";
    const double tolerance = fragment_mass_tolerance_;
    const SpectrumGenerators_ generators = spectrumGenerators_();
    const bool rich = ion_prior_model_ == FragmentIonLikelihoodModel::ContextSet::RICH;
    // Fragment charges of the theoretical ions: b/y ions of a z+ precursor carry at most z - 1 charges.
    const int max_fragment_charge = ion_prior_max_fragment_charge_;
    const auto fragment_charges = [max_fragment_charge](int precursor_charge) {
      return std::max(1, std::min(max_fragment_charge, precursor_charge - 1));
    };
    // The ions of @p sequence (charge 1 to fragment_charges(charge)) matched against peak list @p scan: from the
    // ion types, ordinals and charges of getIonSeriesSpectrum() if the generator supports it (no ion names to write and
    // parse back: the same ions), else from getSpectrum() with the names.
    struct IonScratch
    {
      PeakSpectrum theo;
      std::vector<uint8_t> series;
      std::vector<UInt> ordinals;
      std::vector<Int> charges;
      std::vector<Residue::ResidueType> types;
    };
    const auto match_ions = [&](const FragmentIonLikelihoodModel& model, const TheoreticalSpectrumGenerator& tsg, Size scan,
                                const AASequence& sequence, int charge, IonScratch& scratch, std::vector<FragmentIonLikelihoodModel::Ion>& ions)
    {
      if (tsg.getIonSeriesSpectrum(scratch.theo, scratch.series, scratch.ordinals, scratch.charges, scratch.types, sequence, 1,
                                   fragment_charges(charge)))
      {
        model.matchIons(evidence, scan, scratch.theo, scratch.types, scratch.ordinals, scratch.charges, sequence, charge, tolerance, ppm, ions);
        return;
      }
      scratch.theo.clear(true);
      tsg.getSpectrum(scratch.theo, sequence, 1, fragment_charges(charge));
      model.matchIons(evidence, scan, scratch.theo, sequence, charge, tolerance, ppm, ions);
    };

    // 1. The spectrum of every PSM (scan_index of postProcessHits_); -1 if it has no peaks to match. In parallel, with
    //    the meta values read by registry index (a lookup by name takes the registry's process-wide lock).
    const Size n_ids = peptide_ids.size();
    std::vector<SignedSize> scans(n_ids, -1);
    const UInt mv_scan_index = MetaInfoInterface::metaRegistry().getIndex("scan_index"); // UInt(-1): never registered
    const UInt mv_target_decoy = MetaInfoInterface::metaRegistry().getIndex(Constants::UserParam::TARGET_DECOY);
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      const PeptideIdentification& pi = peptide_ids[i];
      if (pi.getHits().empty() || !pi.metaValueExists(mv_scan_index)) continue;
      const SignedSize scan = static_cast<int>(pi.getMetaValue(mv_scan_index));
      if (scan >= 0 && static_cast<Size>(scan) < spectra.size() && static_cast<Size>(scan) < evidence.size()
          && evidence.peaks(static_cast<Size>(scan)) > 0)
      {
        scans[i] = scan;
      }
    }
    // Best hit by the identification's score orientation (postProcessHits_ sorts the hits, but the
    // annotation must not depend on it).
    const auto best_hit = [](const PeptideIdentification& pi) -> const PeptideHit& {
      const std::vector<PeptideHit>& hits = pi.getHits();
      const bool higher_better = pi.isHigherScoreBetter();
      Size best = 0;
      for (Size h = 1; h < hits.size(); ++h)
      {
        if (higher_better ? hits[h].getScore() > hits[best].getScore() : hits[h].getScore() < hits[best].getScore()) best = h;
      }
      return hits[best];
    };

    // 2. Training PSMs of each half (fold = scan parity): its best target hits at TDC q <= threshold over the
    //    native score, computed within the half, so the other half's labels are not used. Decoys win exact
    //    ties, as in the native yield evaluation, which keeps the selection conservative. The search has decoys
    //    (target-only searches never get here), so (D + 1) / T is an estimate even in a half without a decoy hit.
    //    The candidate rows of the PSMs are made in parallel (label: -1 none, 0 target, 1 decoy) and collected in PSM
    //    order; each half is sorted by a total order (the index decides the last ties), so the order does not depend
    //    on the sort, and the two halves are sorted at the same time.
    struct Row { double score; bool decoy; Size index; };
    std::vector<Row> candidate(n_ids);
    std::vector<signed char> label(n_ids, -1);
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      if (scans[i] < 0) continue;
      const PeptideHit& top = best_hit(peptide_ids[i]);
      if (!top.metaValueExists(mv_target_decoy)) continue;
      const bool decoy = top.getMetaValue(mv_target_decoy).toString() == "decoy";
      candidate[i] = {peptide_ids[i].isHigherScoreBetter() ? top.getScore() : -top.getScore(), decoy, static_cast<Size>(i)};
      label[i] = decoy ? 1 : 0;
    }
    std::array<std::vector<Row>, 2> rows;
    for (Size i = 0; i < n_ids; ++i)
    {
      if (label[i] >= 0) rows[static_cast<Size>(scans[i]) % 2].push_back(candidate[i]);
    }
    std::vector<Row>().swap(candidate);
    std::array<std::vector<Size>, 2> training;
    // (two threads, if this thread has them; declared without OpenMP, too: -fopenmp-simd parses the clause)
    [[maybe_unused]] int sort_threads = 1;
#ifdef _OPENMP
    sort_threads = std::min(2, omp_get_max_threads());
#endif
#pragma omp parallel for schedule(static, 1) num_threads(sort_threads)
    for (int fold = 0; fold < 2; ++fold)
    {
      std::vector<Row>& fold_rows = rows[fold];
      std::sort(fold_rows.begin(), fold_rows.end(), [](const Row& a, const Row& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.decoy != b.decoy) return a.decoy;
        return a.index < b.index;
      });
    }
    for (Size fold = 0; fold < 2; ++fold)
    {
      std::vector<Row>& fold_rows = rows[fold];
      std::vector<double> q(fold_rows.size(), 1.0);
      Size targets = 0, decoys = 0;
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (fold_rows[r].decoy) { ++decoys; }
        else { ++targets; }
        q[r] = static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1));
      }
      for (Size r = fold_rows.size(); r-- > 1;) { q[r - 1] = std::min(q[r - 1], q[r]); }
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (!fold_rows[r].decoy && q[r] <= ion_prior_train_fdr_) training[fold].push_back(fold_rows[r].index);
      }
    }
    const bool trained = training[0].size() >= ion_prior_min_psms_ && training[1].size() >= ion_prior_min_psms_
                         && std::isfinite(tolerance) && tolerance > 0.0;

    // 3. One model per half; its signal: the ions of the training PSMs, its noise: those of their reversed sequences
    //    on the same spectra. Ions are matched in parallel; the integer counts are added up in training order.
    std::array<FragmentIonLikelihoodModel, 2> models{FragmentIonLikelihoodModel(ion_prior_model_),
                                                     FragmentIonLikelihoodModel(ion_prior_model_)};
    if (trained)
    {
      for (Size fold = 0; fold < 2; ++fold)
      {
        const std::vector<Size>& psms = training[fold];
        FragmentIonLikelihoodModel& model = models[fold];
        std::vector<std::vector<FragmentIonLikelihoodModel::Ion>> signal(psms.size()), noise(psms.size());
#pragma omp parallel for schedule(dynamic, 16)
        for (SignedSize t = 0; t < static_cast<SignedSize>(psms.size()); ++t)
        {
          const Size scan = static_cast<Size>(scans[psms[t]]);
          const PeptideHit& hit = best_hit(peptide_ids[psms[t]]);
          const int charge = hit.getCharge();
          const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spectra[scan]);
          IonScratch scratch;
          match_ions(model, tsg, scan, hit.getSequence(), charge, scratch, signal[t]);
          const AASequence reversed = reversedNoiseSequence_(hit.getSequence());
          if (reversed == hit.getSequence()) continue; // no noise hypothesis (a palindrome)
          match_ions(model, tsg, scan, reversed, charge, scratch, noise[t]);
        }
        for (Size t = 0; t < psms.size(); ++t)
        {
          model.addObservations(signal[t], false);
          if (!noise[t].empty()) model.addObservations(noise[t], true);
        }
        model.finalize();
      }
    }
    else
    {
      OPENMS_LOG_WARN << "[ProSE] Ion priors: " << training[0].size() << " and " << training[1].size()
                      << " confident training PSMs in the two halves of the spectra, fewer than " << ion_prior_min_psms_
                      << "; the ion_prior_* features are 0 for this file." << std::endl;
    }
    sw_train.stop();

    // 4. Annotate every hit, scoring it with the model of the other half (cross-fitting). Untrained runs get zeros
    //    so the feature columns stay complete. The names of annotate:ion_prior_features are registered here, always in
    //    this order, before the parallel loop.
    sw_annotate.start();
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const bool write_llr = ion_prior_feature_llr_, write_explained = ion_prior_feature_explained_, write_topk = ion_prior_feature_topk_;
    const UInt mv_llr = write_llr ? registry.registerName(Constants::UserParam::ION_PRIOR_LLR) : 0;
    const UInt mv_explained = write_explained ? registry.registerName(Constants::UserParam::ION_PRIOR_EXPLAINED) : 0;
    const UInt mv_topk = write_topk ? registry.registerName(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED) : 0;
    Size annotated = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : annotated)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      PeptideIdentification& pi = peptide_ids[i];
      const SignedSize scan = scans[i];
      const bool scorable = trained && scan >= 0;
      const FragmentIonLikelihoodModel& model = models[scorable ? 1 - static_cast<Size>(scan) % 2 : 0];
      IonScratch scratch;
      std::vector<FragmentIonLikelihoodModel::Ion> ions;
      for (PeptideHit& hit : pi.getHits())
      {
        FragmentIonLikelihoodModel::Features features;
        if (scorable)
        {
          const int charge = hit.getCharge();
          match_ions(model, generators.forSpectrum(spectra[scan]), static_cast<Size>(scan), hit.getSequence(), charge, scratch, ions);
          features = model.score(ions);
        }
        if (write_llr) hit.setMetaValue(mv_llr, features.log_likelihood_ratio);
        if (write_explained) hit.setMetaValue(mv_explained, features.explained_presence);
        if (write_topk) hit.setMetaValue(mv_topk, features.top_predicted_observed);
        ++annotated;
      }
    }
    sw_annotate.stop();
    OPENMS_LOG_INFO << "[ProSE] Ion priors (" << (rich ? "rich" : "basic") << " model, "
                    << (ion_prior_scored_peaks_ ? "scored" : "all") << " peaks): "
                    << (trained ? "trained on " : "not trained; ") << training[0].size() << " + " << training[1].size()
                    << " PSMs at q <= " << ion_prior_train_fdr_ << " (halves by scan parity, each scored by the other's model) in "
                    << sw_train.getClockTime() << " s; " << annotated << " hits annotated in " << sw_annotate.getClockTime()
                    << " s; peak lists " << evidence.totalPeaks() << " peaks of " << evidence.size() << " spectra, "
                    << static_cast<double>(evidence.memoryUsage()) / (1024.0 * 1024.0) << " MiB." << std::endl;

    // 5. Percolator sees the new columns; the search parameters record every setting that determines the features
    //    (model, peaks, fragment charges, training threshold and minimum) and the training set.
    if (protein_ids.empty()) return;
    ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
    StringList features;
    if (search_parameters.metaValueExists("extra_features"))
    {
      const std::string existing = search_parameters.getMetaValue("extra_features").toString();
      if (!existing.empty()) features = ListUtils::create<std::string>(existing);
    }
    StringList written;
    if (write_llr) written.push_back(Constants::UserParam::ION_PRIOR_LLR);
    if (write_explained) written.push_back(Constants::UserParam::ION_PRIOR_EXPLAINED);
    if (write_topk) written.push_back(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED);
    for (const std::string& name : written)
    {
      if (std::find(features.begin(), features.end(), name) == features.end()) features.push_back(name);
    }
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
    search_parameters.setMetaValue("annotate:self_trained_ion_priors", "true");
    search_parameters.setMetaValue("annotate:ion_prior_features", ListUtils::concatenate(written, ","));
    search_parameters.setMetaValue("annotate:ion_prior_model", rich ? "rich" : "basic");
    search_parameters.setMetaValue("annotate:ion_prior_peaks", ion_prior_scored_peaks_ ? "scored" : "all");
    search_parameters.setMetaValue("annotate:ion_prior_max_fragment_charge", ion_prior_max_fragment_charge_);
    search_parameters.setMetaValue("annotate:ion_prior_train_fdr", ion_prior_train_fdr_);
    search_parameters.setMetaValue("annotate:ion_prior_min_psms", static_cast<int>(ion_prior_min_psms_));
    search_parameters.setMetaValue("ion_prior:trained", trained ? "true" : "false");
    search_parameters.setMetaValue("ion_prior:training_psms", static_cast<int>(training[0].size() + training[1].size()));
    search_parameters.setMetaValue("ion_prior:fold_training_psms",
                                   IntList{static_cast<int>(training[0].size()), static_cast<int>(training[1].size())});
  }
} // namespace OpenMS
