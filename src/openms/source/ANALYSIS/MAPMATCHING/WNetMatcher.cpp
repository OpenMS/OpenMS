// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/WNetMatcher.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <wnetalign/aligner.hpp>
#include <wnetalign/spectrum.hpp>

using namespace std;

namespace OpenMS
{

  WNetMatcher::DistanceMetric WNetMatcher::metricFromString(const string& s)
  {
    if (s == "L1") return DistanceMetric::L1;
    if (s == "L2") return DistanceMetric::L2;
    if (s == "LINF") return DistanceMetric::LINF;
    OPENMS_LOG_WARN << "WNetMatcher::metricFromString: unrecognized metric '" << s
                    << "', falling back to LINF.\n";
    return DistanceMetric::LINF;
  }

  namespace
  {
    ::DistanceMetric toWNetMetric(WNetMatcher::DistanceMetric m)
    {
      switch (m)
      {
        case WNetMatcher::DistanceMetric::L1:   return ::DistanceMetric::L1;
        case WNetMatcher::DistanceMetric::L2:   return ::DistanceMetric::L2;
        case WNetMatcher::DistanceMetric::LINF: return ::DistanceMetric::LINF;
        default:                                return ::DistanceMetric::LINF;
      }
    }
  } // anonymous namespace

  WNetMatcher::MatchResult WNetMatcher::match(
    const vector<array<double, 2>>& positions_a,
    const vector<double>& intensities_a,
    const vector<array<double, 2>>& positions_b,
    const vector<double>& intensities_b,
    DistanceMetric metric,
    double max_distance,
    double trash_cost)
  {
    MatchResult result;
    result.cost = 0.0;

    if (positions_a.size() != intensities_a.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "positions_a and intensities_a must have the same length");
    }
    if (positions_b.size() != intensities_b.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "positions_b and intensities_b must have the same length");
    }

    if (positions_a.empty() || positions_b.empty())
    {
      return result;
    }

    // Build spectra: spec_a is the empirical set, spec_b the theoretical target.
    Spectrum<2> spec_a(positions_a, intensities_a);
    Spectrum<2> spec_b(positions_b, intensities_b);

    // Construct aligner for empirical (spec_a) vs. one theoretical target (spec_b).
    vector<Spectrum<2>*> theoretical = {&spec_b};
    WNetAligner<2> aligner(spec_a, theoretical, toWNetMetric(metric), max_distance, trash_cost);
    // set_point({1.0}) fixes a single unit-weight consensus and triggers the network-flow solve.
    aligner.set_point({1.0});

    // consensus_for_target(0) returns paired index lists: emp_ids[k] in spec_a matched to theo_ids[k] in spec_b.
    auto [emp_ids, theo_ids] = aligner.consensus_for_target(0);

    // total_cost() reads the optimal transport cost computed by set_point().
    result.cost = aligner.total_cost();
    result.matched_pairs.reserve(emp_ids.size());
    for (Size i = 0; i < emp_ids.size(); ++i)
    {
      result.matched_pairs.emplace_back(
        static_cast<Size>(emp_ids[i]),
        static_cast<Size>(theo_ids[i]));
    }
    return result;
  }

} // namespace OpenMS
