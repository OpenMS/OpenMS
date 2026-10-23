// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXDatasetExport.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <deque>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <unordered_map>

namespace OpenMS
{
  namespace
  {
    /// Whether a higher value of the score @p name (a QPXDatasetExport::scoreName()) is better; none if not known
    std::optional<bool> higherBetter(const std::string& name)
    {
      if (name == "global_qvalue" || name == "percolator_qvalue" || name == "percolator_pep") { return false; }
      if (name == "percolator_score" || name == "ln_hyperscore") { return true; }
      return std::nullopt;
    }
  }

  std::string QPXDatasetExport::scoreName(const std::string& name)
  {
    static const std::map<std::string, std::string> known = {
      {"q-value", "global_qvalue"},
      {"MS:1001491", "percolator_qvalue"}, {"percolator_q_value", "percolator_qvalue"},
      {"MS:1001492", "percolator_score"}, {"svm", "percolator_score"},
      {"MS:1001493", "percolator_pep"},
      {"ln(hyperscore)", "ln_hyperscore"}, {"ln(hyperscore)_score", "ln_hyperscore"},
      {"Posterior Error Probability", ""}, {"Posterior Error Probability_score", ""},
      {"IsotopeError", "isotope_error"}};
    const auto found = known.find(name);
    if (found != known.end()) { return found->second; }
    std::string out;
    for (Size i = 0; i < name.size(); ++i)
    {
      const unsigned char c = static_cast<unsigned char>(name[i]);
      const unsigned char before = i > 0 ? static_cast<unsigned char>(name[i - 1]) : 0;
      if (std::isupper(c) && (std::islower(before) || std::isdigit(before))) { out += '_'; }
      if (std::isalnum(c)) { out += static_cast<char>(std::tolower(c)); }
      else if (!out.empty() && out.back() != '_') { out += '_'; }
    }
    while (!out.empty() && out.back() == '_') { out.pop_back(); }
    return out;
  }

  std::vector<QPXDatasetExport::Column> QPXDatasetExport::psmColumns(const PeptideIdentificationList& peptide_ids,
                                                                     const std::vector<bool>& has_spectrum,
                                                                     bool isotope_observed_minus_theoretical)
  {
    if (!has_spectrum.empty() && has_spectrum.size() != peptide_ids.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "psmColumns: has_spectrum has " + std::to_string(has_spectrum.size()) + " entries for "
        + std::to_string(peptide_ids.size()) + " peptide identifications");
    }

    auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
    Column scores = column(QPXPSMSchema::ADDITIONAL_SCORES, Column::Type::SCORE_LIST);
    Column reference = column("spectrum_reference", Column::Type::UTF8);
    Column isotope = column("isotope_error", Column::Type::INT16);
    Column mappings = column("protein_mappings", Column::Type::STRUCT_LIST);
    mappings.struct_fields = {{"accession", Column::Type::UTF8, false}, {"start", Column::Type::INT32}, {"end", Column::Type::INT32},
                              {"aa_before", Column::Type::UTF8}, {"aa_after", Column::Type::UTF8}};

    // The isotope error as QPX has it (bigbio/qpx#355: the 13C peaks from the monoisotopic to the selected peak, i.e.
    // observed minus theoretical), in the column and in 'additional_scores'
    const int isotope_sign = isotope_observed_minus_theoretical ? 1 : -1;
    // meta values that are no scores (columns of their own, or attributes of the spectrum)
    static const std::set<std::string> not_scores = {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
                                                     "FAIMS_CV", "scan", "reference_file_name", "id_merge_index"};
    auto residue = [](char aa) -> std::optional<std::string>
    {
      if (aa == PeptideEvidence::N_TERMINAL_AA || aa == PeptideEvidence::C_TERMINAL_AA) { return std::string("-"); }
      if (aa == PeptideEvidence::UNKNOWN_AA) { return std::nullopt; }
      return std::string(1, aa);
    };
    auto position = [](Int p) -> std::optional<Int32>
    {
      if (p == PeptideEvidence::UNKNOWN_POSITION || p < 0) { return std::nullopt; }
      return static_cast<Int32>(p + 1); // 1-based
    };

    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      const PeptideIdentification& id = peptide_ids[i];
      if (id.getHits().empty()) { continue; } // no QPX row
      const PeptideHit& hit = id.getHits()[0];

      std::vector<Column::Score> row_scores;
      std::set<std::string> seen;
      auto add = [&](const std::string& key, double value, std::optional<bool> higher_better)
      {
        const std::string name = scoreName(key);
        if (name.empty() || !seen.insert(name).second) { return; }
        const std::optional<bool> known = higherBetter(name);
        if (name == "isotope_error" && value != 0.0) { value *= isotope_sign; } // (no -0)
        row_scores.push_back({name, value, known ? known : higher_better});
      };
      add(id.getScoreType(), hit.getScore(), id.isHigherScoreBetter());
      std::vector<std::string> keys;
      hit.getKeys(keys);
      for (const std::string& key : keys)
      {
        if (not_scores.count(key) > 0) { continue; }
        const DataValue& value = hit.getMetaValue(key);
        if (value.valueType() == DataValue::INT_VALUE || value.valueType() == DataValue::DOUBLE_VALUE)
        {
          add(key, static_cast<double>(value), std::nullopt);
        }
      }
      scores.score_lists.push_back(std::move(row_scores));

      // (null: the mz view has no row of its spectrum, or no native ID)
      const bool spectrum = has_spectrum.empty() || has_spectrum[i];
      reference.utf8_values.push_back(spectrum ? id.getSpectrumReference() : std::string());
      reference.null_rows.push_back(!spectrum || id.getSpectrumReference().empty());

      const std::string isotope_key = hit.metaValueExists("isotope_error") ? "isotope_error" : "IsotopeError";
      isotope.int64_values.push_back(hit.metaValueExists(isotope_key) ? isotope_sign * static_cast<int>(hit.getMetaValue(isotope_key)) : 0);
      isotope.null_rows.push_back(!hit.metaValueExists(isotope_key));

      std::vector<Column::StructEntry> row_mappings;
      for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
      {
        row_mappings.push_back({evidence.getProteinAccession(), position(evidence.getStart()), position(evidence.getEnd()),
                                residue(evidence.getAABefore()), residue(evidence.getAAAfter())});
      }
      mappings.struct_lists.push_back(std::move(row_mappings));
    }
    return {std::move(scores), std::move(reference), std::move(isotope), std::move(mappings)};
  }

  std::vector<Size> QPXDatasetExport::keptPeaks(const MSSpectrum& spectrum, const std::vector<PeptideHit::PeakAnnotation>& annotations,
                                                int precursor_charge, const PeakFilter& filter)
  {
    std::vector<bool> keep(spectrum.size(), filter.top == 0);
    auto indexOf = [&spectrum](double mz) -> Size // the peak nearest to mz
    {
      const auto above = std::lower_bound(spectrum.begin(), spectrum.end(), mz,
                                          [](const Peak1D& p, double m) { return p.getMZ() < m; });
      Size i = static_cast<Size>(above - spectrum.begin());
      if (i == spectrum.size() || (i > 0 && mz - spectrum[i - 1].getMZ() < above->getMZ() - mz)) { --i; }
      return i;
    };
    if (filter.top > 0 && !spectrum.empty())
    {
      // Jumping windows as WindowMower's (a window starts at its first peak; the next at the first peak past it),
      // but the last window keeps the full quota too: WindowMower scales it by the window's width, which would
      // drop the highest-m/z peaks of nearly every spectrum.
      std::vector<Size> order;
      for (Size begin = 0, end = 0; begin < spectrum.size(); begin = end)
      {
        while (end < spectrum.size() && spectrum[end].getMZ() - spectrum[begin].getMZ() < filter.window) { ++end; }
        order.resize(end - begin);
        std::iota(order.begin(), order.end(), begin);
        const Size n = std::min(filter.top, order.size());
        std::partial_sort(order.begin(), order.begin() + n, order.end(), [&spectrum](Size a, Size b)
        {
          return spectrum[a].getIntensity() > spectrum[b].getIntensity() || (spectrum[a].getIntensity() == spectrum[b].getIntensity() && a < b);
        });
        for (Size k = 0; k < n; ++k) { keep[order[k]] = true; }
      }
      int max_charge = filter.deisotoped ? filter.max_charge : 1;
      if (filter.deisotoped && filter.cap_at_precursor && precursor_charge > 0) { max_charge = std::min(max_charge, precursor_charge); }
      for (const PeptideHit::PeakAnnotation& a : annotations)
      {
        for (int z = 1; z <= std::max(1, max_charge); ++z)
        {
          const double mz = (a.mz + (z - 1) * Constants::PROTON_MASS_U) / z;
          const Size i = indexOf(mz);
          if (std::fabs(spectrum[i].getMZ() - mz) <= 0.01) { keep[i] = true; }
        }
      }
    }
    std::vector<Size> kept;
    for (Size i = 0; i < keep.size(); ++i) { if (keep[i]) { kept.push_back(i); } }
    return kept;
  }

  QPXDatasetExport::MzView QPXDatasetExport::mzView(const MSExperiment& spectra, const PeptideIdentificationList& peptide_ids,
                                                    const std::string& run_file_name, const PeakFilter& filter)
  {
    MzView result;
    result.has_spectrum.assign(peptide_ids.size(), false);

    // -- the MS2 spectra of the PSMs (several, if a native ID is not unique), sorted by m/z --
    std::unordered_map<std::string, std::vector<Size>> psms_of_spectrum; // native ID -> PSMs
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (!peptide_ids[i].getHits().empty()) { psms_of_spectrum[peptide_ids[i].getSpectrumReference()].push_back(i); }
    }
    std::deque<MSSpectrum> sorted; // the sorted copies of unsorted spectra (stable addresses)
    std::vector<const MSSpectrum*> of_psms;
    for (const MSSpectrum& spectrum : spectra)
    {
      if (spectrum.getMSLevel() != 2 || psms_of_spectrum.count(spectrum.getNativeID()) == 0) { continue; }
      if (spectrum.isSorted()) { of_psms.push_back(&spectrum); continue; }
      sorted.push_back(spectrum);
      sorted.back().sortByPosition();
      of_psms.push_back(&sorted.back());
    }
    std::unordered_map<std::string, std::vector<const MSSpectrum*>> spectra_of_id;
    for (const MSSpectrum* spectrum : of_psms) { spectra_of_id[spectrum->getNativeID()].push_back(spectrum); }
    // the spectrum of PSM i: by its native ID, of several the one at its RT; null if none
    auto spectrumOf = [&](Size i) -> const MSSpectrum*
    {
      const auto found = spectra_of_id.find(peptide_ids[i].getSpectrumReference());
      if (found == spectra_of_id.end()) { return nullptr; }
      if (found->second.size() == 1) { return found->second.front(); }
      for (const MSSpectrum* spectrum : found->second)
      {
        if (std::fabs(spectrum->getRT() - peptide_ids[i].getRT()) <= 1e-6) { return spectrum; }
      }
      return nullptr;
    };

    // -- the mz view: a row per spectrum, in the order of its first PSM --
    using Value = Column::StructValue;
    auto column = [](const std::string& name, Column::Type type, bool nullable = true)
    {
      Column c; c.name = name; c.type = type; c.nullable = nullable; return c;
    };
    Column id = column("id", Column::Type::UTF8, false);
    Column run_file = column("run_file_name", Column::Type::UTF8);
    Column scan = column("scan", Column::Type::INT32_LIST);
    Column ms_level = column("ms_level", Column::Type::INT32, false);
    Column centroid = column("centroid", Column::Type::BOOL, false);
    Column start_time = column("scan_start_time", Column::Type::FLOAT, false);
    Column mobility = column("inverse_ion_mobility", Column::Type::FLOAT);
    Column injection = column("ion_injection_time", Column::Type::FLOAT);
    Column tic = column("total_ion_current", Column::Type::FLOAT, false);
    Column precursors = column("precursors", Column::Type::STRUCT_LIST);
    precursors.struct_fields = {{"selected_ion_mz", Column::Type::FLOAT, false}, {"selected_ion_charge", Column::Type::INT32},
                                {"selected_ion_intensity", Column::Type::FLOAT}, {"isolation_window_target", Column::Type::FLOAT},
                                {"isolation_window_lower", Column::Type::FLOAT}, {"isolation_window_upper", Column::Type::FLOAT},
                                {"spectrum_ref", Column::Type::UTF8}};
    Column mz = column("mz", Column::Type::FLOAT_LIST, false);
    Column intensity = column("intensity", Column::Type::FLOAT_LIST, false);
    Column cv_params = column("cv_params", Column::Type::STRUCT_LIST);
    cv_params.struct_fields = {{"cv_name", Column::Type::UTF8, false}, {"cv_value", Column::Type::UTF8, false}};
    auto optional_float = [](Column& c, double v, bool valid)
    {
      c.double_values.push_back(valid ? v : 0.0);
      c.null_rows.push_back(!valid);
    };
    // the spectrum with a row of each native ID: one row per native ID ('id' is the mz view's primary key), of
    // several spectra of one native ID (at other RTs) the first's; the PSMs of the others have none
    std::map<std::string, const MSSpectrum*> written;
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (peptide_ids[i].getHits().empty()) { continue; }
      const MSSpectrum* spectrum = spectrumOf(i);
      if (spectrum == nullptr) { continue; }
      if (!written.emplace(spectrum->getNativeID(), spectrum).second) { continue; } // its native ID has a row
      // its peaks: those the filter keeps for any of its PSMs (each keeps its own annotated peaks), with the
      // deisotoper's charge cap of the search's preprocessing: the spectrum's recorded precursor charge (0: unknown,
      // then up to filter.max_charge), not the charge of a PSM (a spectrum of unknown charge was deisotoped up to
      // that maximum, whatever the charge of its PSMs)
      const int precursor_charge = spectrum->getPrecursors().empty() ? 0 : spectrum->getPrecursors()[0].getCharge();
      std::vector<bool> keep(spectrum->size(), false);
      for (Size j : psms_of_spectrum[peptide_ids[i].getSpectrumReference()])
      {
        if (spectrumOf(j) != spectrum) { continue; }
        const PeptideHit& hit = peptide_ids[j].getHits()[0];
        for (Size k : keptPeaks(*spectrum, hit.getPeakAnnotations(), precursor_charge, filter)) { keep[k] = true; }
      }
      std::vector<float> peak_mz, peak_intensity;
      for (Size k = 0; k < keep.size(); ++k)
      {
        if (!keep[k]) { continue; }
        peak_mz.push_back(static_cast<float>((*spectrum)[k].getMZ()));
        peak_intensity.push_back((*spectrum)[k].getIntensity());
      }
      result.peaks_read += spectrum->size();
      result.peaks_kept += peak_mz.size();

      id.utf8_values.push_back(spectrum->getNativeID());
      run_file.utf8_values.push_back(run_file_name);
      scan.int32_lists.push_back(ArrowIOHelpers::qpxScanComponents(spectrum->getNativeID()));
      ms_level.int64_values.push_back(static_cast<Int64>(spectrum->getMSLevel()));
      SpectrumSettings::SpectrumType type = spectrum->getType(false);
      if (type == SpectrumSettings::SpectrumType::UNKNOWN) { type = spectrum->getType(true); }
      centroid.int64_values.push_back(type == SpectrumSettings::SpectrumType::CENTROID ? 1 : 0);
      start_time.double_values.push_back(spectrum->getRT() / 60.0);
      optional_float(mobility, spectrum->getDriftTime(),
                     spectrum->getDriftTimeUnit() == DriftTimeUnit::VSSC && spectrum->getDriftTime() >= 0.0);
      const bool has_injection = !spectrum->getAcquisitionInfo().empty()
                                 && spectrum->getAcquisitionInfo()[0].metaValueExists("MS:1000927");
      optional_float(injection, has_injection ? static_cast<double>(spectrum->getAcquisitionInfo()[0].getMetaValue("MS:1000927")) : 0.0,
                     has_injection);
      tic.double_values.push_back(spectrum->calculateTIC());
      std::vector<Column::StructEntry> ions;
      for (const Precursor& p : spectrum->getPrecursors())
      {
        ions.push_back({p.getMZ(), p.getCharge() > 0 ? Value(static_cast<Int32>(p.getCharge())) : Value(),
                        p.getIntensity() > 0 ? Value(static_cast<double>(p.getIntensity())) : Value(), Value(),
                        p.getIsolationWindowLowerOffset() > 0 ? Value(p.getIsolationWindowLowerOffset()) : Value(),
                        p.getIsolationWindowUpperOffset() > 0 ? Value(p.getIsolationWindowUpperOffset()) : Value(), Value()});
      }
      precursors.struct_lists.push_back(std::move(ions));
      mz.float_lists.push_back(std::move(peak_mz));
      intensity.float_lists.push_back(std::move(peak_intensity));
      cv_params.struct_lists.emplace_back();
      cv_params.null_rows.push_back(true);
    }
    const Size mz_rows = id.utf8_values.size();
    result.table = ArrowIOHelpers::tableOfColumns({id, run_file, scan, ms_level, centroid, start_time, mobility, injection, tic,
                                                   precursors, mz, intensity, cv_params}, mz_rows);

    // Each PSM against the rows, once all are made (a later PSM may make the row of a native ID): a PSM has its
    // spectrum's row only if the row of its native ID is that spectrum's. Otherwise (no spectrum of its native ID and
    // RT, no table, or another spectrum of its native ID has the row) it has none: a join on its native ID would give
    // no peaks or another spectrum's
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      if (peptide_ids[i].getHits().empty()) { continue; }
      const auto row = written.find(peptide_ids[i].getSpectrumReference());
      if (row == written.end()) { ++result.without_spectrum; } // no spectrum of its native ID (and RT), and no row to join
      else if (row->second != spectrumOf(i)) { ++result.unresolved; }
      else { result.has_spectrum[i] = result.table != nullptr; }
    }
    return result;
  }
} // namespace OpenMS
