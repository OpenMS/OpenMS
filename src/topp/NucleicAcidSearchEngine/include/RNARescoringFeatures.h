// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <optional>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Percolator feature rows for the oligonucleotide-spectrum matches of NASE.

  Every row has the same columns (featureNames()); which columns exist depends
  only on the Options, never on a hit, so targets and decoys are described
  identically. Rows are computed for all hits of one spectrum together, because
  the delta features compare each hit with its competitors (computeRow() for one hit):

  - native_score: the search score the hits were ranked by
  - delta_next: score minus the best lower score of a hit with a different
    sequence (minus 0 if there is none; native scores are positive)
  - delta_best: score minus the best score of the spectrum (0 for the best hit)
  - abs_ppm_error, isotope_offset, charge_<z> (one-hot over |z| in
    [min_charge, max_charge]), length, modifications, adduct
  - frac_<series>: distinct matched ordinals of one ion series / (length - 1)
  - matched_ions: distinct matched (series, ordinal) pairs / (series * (length - 1))
  - matched_intensity: intensity of the distinct matched peaks / total spectrum intensity
  - longest_series: longest run of consecutive ordinals of any series / (length - 1)
  - coverage: internal backbone positions explained by any series / (length - 1)
  - route_nodes / route_intensity_ratio / route_prefix_length / route_suffix_length
    (Options::route: the route score is the search score, so native_score),
    graph_tag_length (Options::graph), diagnostic_supported_mods (Options::diagnostics:
    modifications permitted by the spectrum's diagnostic ions; the unsupported ones
    are modifications - diagnostic_supported_mods),
    flash_corrected / flash_cosine / flash_snr (Options::flash)

  Positions follow NASE's coverage rule: prefix ion n (a, b, c, d, a-B)
  explains the bond after residue n, suffix ion n (w, x, y, z) the bond after
  residue length - n. Annotations whose ion series is not configured, or which
  cannot be parsed, are ignored.
*/
class RNARescoringFeatures
{
public:
  struct Options
  {
    /// Ion series of the theoretical spectra (NASE 'fragment:ions'); one frac_<series> column each, in this order.
    StringList ion_series;
    /// Range of absolute precursor charges of the charge_<z> one-hot columns.
    Size min_charge = 1, max_charge = 1;
    bool route = false;
    bool graph = false;
    bool diagnostics = false;
    bool flash = false;
  };

  /// Per-spectrum context shared by all hits of a spectrum.
  struct SpectrumInput
  {
    /// Total intensity of the spectrum the hits were scored against (same scale as the annotation intensities).
    double total_intensity = 0.0;
    /// Modification codes permitted by diagnostic evidence (Options::diagnostics).
    std::set<std::string> diagnostic_allowed;
    /// FLASH precursor cosine and SNR, if the precursor was corrected (Options::flash).
    std::optional<double> flash_cosine, flash_snr;
  };

  /// One hit. Hits of a spectrum may be passed in any order.
  struct HitInput
  {
    double native_score = 0.0;
    NASequence sequence;
    double precursor_error_ppm = 0.0;
    Size isotope_offset = 0;
    Int charge = 0;
    bool adduct = false;
    const std::vector<PeptideHit::PeakAnnotation>* annotations = nullptr;
    // Options::route
    Size route_nodes = 0;
    double route_matched_intensity = 0.0;
    double route_local_intensity = 0.0;
    Size route_prefix_length = 0;
    Size route_suffix_length = 0;
    // Options::graph
    Size graph_tag_length = 0;
  };

  explicit RNARescoringFeatures(Options options);

  /// Column names, aligned with every row of computeSpectrum().
  const StringList& featureNames() const { return names_; }

  /**
    @brief The feature row of hit @p index of one spectrum; the other @p hits only give its score margins.

    Linear in the number of hits. Equal to computeSpectrum(hits, spectrum)[index].
  */
  std::vector<double> computeRow(const std::vector<HitInput>& hits, Size index, const SpectrumInput& spectrum) const;

  /// One feature row per hit of one spectrum, in the order of @p hits (quadratic in the number of hits; use
  /// computeRow() if only some rows are needed).
  std::vector<std::vector<double>> computeSpectrum(const std::vector<HitInput>& hits,
                                                   const SpectrumInput& spectrum) const;

  /// Ion series and ordinal of a NASE fragment annotation ("a2-B" -> ("a-B", 2), "w3" -> ("w", 3)).
  static std::optional<std::pair<std::string, Size>> parseAnnotation(const std::string& annotation);

  /// Whether @p series is a 5' (prefix) ion series.
  static bool isPrefixSeries(const std::string& series);

private:
  /// The row of @p hit, given the best lower score of a different sequence (@p next) and the spectrum's best score
  std::vector<double> buildRow_(const HitInput& hit, double next, double best, const SpectrumInput& spectrum) const;

  void appendIonFeatures_(const HitInput& hit, const SpectrumInput& spectrum, std::vector<double>& row) const;

  Options options_;
  StringList names_;
};
} // namespace OpenMS
