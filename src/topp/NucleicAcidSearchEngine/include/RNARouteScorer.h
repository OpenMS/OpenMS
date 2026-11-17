// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include "RNATagCandidateFilter.h"

#include <OpenMS/METADATA/PeptideHit.h>
#include <map>

namespace OpenMS
{
/**
  @brief Rank candidates using terminally anchored prefix and suffix ladders.

  Routes start at ion 1 or 2, retain one series and orientation, and bridge at
  most one missing rung. Charge hypotheses share observed-peak identities.
  The longest route per series is selected, with weighted intensity breaking
  ties. The best pair of these prefix/suffix routes maximizes
  n*(1+I_matched/I_local), crediting an observed peak once. Gaussian mass
  weighting follows NASE's existing scorer.
  Local background is the union of the strongest configured peaks within a
  neutral-mass radius of each matched node. These explicit choices are NASE's
  implementation; the paper does not publish the corresponding numerical rules.
*/
class RNARouteScorer
{
public:
  struct Options
  {
    double tolerance = 20.0;
    bool tolerance_ppm = true;
    /// Gaussian weighting scale in fragment-tolerance units; zero inherits tolerance.
    double mass_accuracy_tolerance = 0.0;
    double local_radius = 50.0;
    Size local_peaks = 10;
    Size max_step = 2;
    Size min_nodes = 2;
    Size max_unknown_charge = 3;
    /// Maximum prefix states constructed across ion series in one match query.
    Size max_states = 200000;
  };
  struct Result
  {
    double score = 0;
    Size nodes = 0;
    double matched_intensity = 0;
    double local_intensity = 0;
    Size prefix_length = 0;
    Size suffix_length = 0;
    std::string prefix_scaffold, suffix_scaffold;
    std::vector<PeptideHit::PeakAnnotation> annotations;
  };

  struct PreparedSpectrum
  {
    /// Per observed peak and allowed charge: unique high-intensity neighbors in neutral-mass space.
    std::vector<std::map<Size, std::vector<Size>>> backgrounds;
  };

  struct MatchedNode
  {
    Size peak, position, charge;
    double weighted, neutral;
    std::string annotation;
  };
  struct Route
  {
    std::vector<MatchedNode> matches;
    Size length = 0;
    double weighted_intensity = 0;
    char series = 0;
    std::string scaffold;
  };
  struct PreparedRoutes
  {
    std::vector<Route> routes;
  };

  /**
    @brief Match terminal ladders once for a reusable route scaffold.
    @param[in] query Annotated scaffold theory.
    @param[in] experimental Read-only processed MS2.
    @param[in] prepared Charge-consistent local contexts for that spectrum.
    @param[in] negative_mode Whether experimental ions are negative.
    @param[in] options Fragment tolerances and route rules.
    @param[in] scaffold Sequence label retained for evidence provenance.
    @return Longest supported route per conventional series.
    @throws Exception::InvalidValue if the prefix-state limit is exceeded; no partial route set is returned.
  */
  static PreparedRoutes match(const RNATagCandidateFilter::PreparedQuery& query,
                              const MSSpectrum& experimental,
                              const PreparedSpectrum& prepared,
                              bool negative_mode,
                              const Options& options,
                              const std::string& scaffold = {});

  /**
    @brief Combine intact and diagnostic-permitted labile-completion route evidence.
    @param[in] evidence Matched route sets; intact evidence takes precedence in exact ties.
    @param[in] experimental Original experimental peak identities.
    @param[in] prepared Cached local backgrounds for that spectrum.
    @param[in] negative_mode Whether experimental ions are negative.
    @return Best prefix/suffix pair after selecting the longest route per series across scaffolds.
  */
  static Result
  score(const std::vector<const PreparedRoutes*>& evidence, const MSSpectrum& experimental, const PreparedSpectrum& prepared, bool negative_mode);

  /**
    @brief Cache charge-consistent local backgrounds once per experimental spectrum.
    @param[in] experimental Processed MS2, sorted by m/z.
    @param[in] negative_mode Whether experimental ions are negative.
    @param[in] options Charge hypotheses and local-background settings.
    @return Reusable contexts with original peak identities.
  */
  static PreparedSpectrum prepareSpectrum(const MSSpectrum& experimental, bool negative_mode, const Options& options);

  /**
    @brief Score with previously prepared experimental contexts.
    @param[in] query Annotated theoretical ions.
    @param[in] experimental The same read-only spectrum used to prepare contexts.
    @param[in] prepared Cached charge hypotheses and local backgrounds.
    @param[in] negative_mode Whether experimental ions are negative.
    @param[in] options Mass tolerances and route settings used during preparation.
    @return Selected nonredundant route evidence, score and annotations.
  */
  static Result score(const RNATagCandidateFilter::PreparedQuery& query,
                      const MSSpectrum& experimental,
                      const PreparedSpectrum& prepared,
                      bool negative_mode,
                      const Options& options);

  /**
    @brief Score a candidate against annotated theoretical ions and read-only MS2.
    @param[in] query Theoretical annotations prepared by RNATagCandidateFilter.
    @param[in] experimental Legacy processed spectrum, sorted by m/z.
    @param[in] negative_mode Whether the spectrum contains negative ions.
    @param[in] options Mass tolerance, route and local-intensity settings.
    @return Selected nonredundant route evidence, score and annotations.
  */
  static Result score(const RNATagCandidateFilter::PreparedQuery& query, const MSSpectrum& experimental, bool negative_mode, const Options& options);
};
} // namespace OpenMS
