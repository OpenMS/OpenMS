// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{

  namespace
  {
    // Helper 1: Charge compatibility check
    bool isChargeCompatible(const ConsensusFeature& a, const ConsensusFeature& b)
    {
      const Int charge_1 = a.getCharge();
      const Int charge_2 = b.getCharge();

      if (charge_1 != 0 && charge_2 != 0 && std::abs(charge_1) != std::abs(charge_2))
      {
        return false;
      }

      const Int effective_charge = (charge_1 != 0) ? std::abs(charge_1) : std::abs(charge_2);
      if (effective_charge == 0)
      {
        return false;
      }

      return true;
    }

    // Helper 2: Isotope shift math and tolerance check
    bool isIsotopeShiftMatch(const ConsensusFeature& a, const ConsensusFeature& b, int max_shift, double c13_mass, double mz_tol, bool mz_ppm)
    {
      const Int charge_1 = a.getCharge();
      const Int charge_2 = b.getCharge();
      const Int effective_charge = (charge_1 != 0) ? std::abs(charge_1) : std::abs(charge_2);

      double mz_diff = std::abs(a.getMZ() - b.getMZ());

      for (int k = 1; k <= max_shift; ++k)
      {
        double expected_shift = (k * c13_mass) / effective_charge;
        double shift_diff = std::abs(mz_diff - expected_shift);
        double tolerance = mz_ppm ? (mz_tol * a.getMZ() / 1e6) : mz_tol;

        if (shift_diff <= tolerance)
        {
          return true;
        }
      }
      return false;
    }

    // Helper 3: Merge, annotate, log, and erase
    void mergeAndAnnotateRescue(ConsensusFeature& target, ConsensusMap::iterator& src_it, int candidate_map_index, ConsensusMap& out)
    {
      // 1. Merge the feature
      target.insert(*src_it->getFeatures().begin());

      // 2. Annotate the rescue for downstream data integrity
      target.setMetaValue("isotope_shift_rescued", "true");

      // 3. Log the successful rescue event
      OPENMS_LOG_DEBUG << "FeatureGroupingAlgorithmQT: Isotope-shift rescue applied. Merged map "
                       << candidate_map_index << " into consensus feature.\n";

      // 4. Remove the source feature from the map
      src_it = out.erase(src_it);
      --src_it;
    }
  }

  FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmQT");
    defaults_.insert("", QTClusterFinder().getParameters());

    defaults_.setValue("enable_isotope_shift_fallback", "false", "If true, performs a second linking pass to rescue features missing a monoisotopic peak.");
    defaults_.setValidStrings("enable_isotope_shift_fallback", {"true","false"});

    defaults_.setValue("max_isotope_shift", 2, "Maximum number of isotopic peaks to shift when matching missing M+0 features (e.g., 1 for M+1, 2 for M+2).");
    defaults_.setMinInt("max_isotope_shift", 1);

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT() = default;

  template <typename MapType>
  void FeatureGroupingAlgorithmQT::group_(const vector<MapType>& maps, ConsensusMap& out)
  {
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "At least two maps must be given!");
    }

    QTClusterFinder cluster_finder;
    cluster_finder.setParameters(param_.copy("", true));

    // PASS 1: Standard exact m/z clustering
    cluster_finder.run(maps, out);

    // PASS 2: Isotopic shift rescue
    if (param_.getValue("enable_isotope_shift_fallback").toBool())
    {
      const Int max_shift = static_cast<Int>(param_.getValue("max_isotope_shift"));
      constexpr double c13_mass = 1.0033548;

      double mz_tol = param_.getValue("distance_MZ:max_difference");
      bool mz_ppm = param_.getValue("distance_MZ:unit").toString() == "ppm";
      double rt_tol = param_.getValue("distance_RT:max_difference");

      // Greedy pass to merge unlinked features (singletons)
      for (auto it1 = out.begin(); it1 != out.end(); ++it1)
      {
        // Only look at features that failed to link in Pass 1
        if (it1->getFeatures().size() != 1) continue;

        for (auto it2 = it1 + 1; it2 != out.end(); ++it2)
        {
          if (it2->getFeatures().size() != 1) continue;

          // Ensure candidate map is not already represented in the merged group
          const auto candidate_map_index = it2->getFeatures().begin()->getMapIndex();
          const bool map_already_present = std::any_of(
            it1->getFeatures().begin(), it1->getFeatures().end(),
            [candidate_map_index](const auto& handle)
            {
              return handle.getMapIndex() == candidate_map_index;
            });

          if (map_already_present)
          {
            continue;
          }

          // Basic RT check using configured tolerance
          if (std::abs(it1->getRT() - it2->getRT()) > rt_tol) continue;

          // Delegate to Helper: Charges must match when both are known
          if (!isChargeCompatible(*it1, *it2)) continue;

          // Delegate to Helper: Check for isotope shifts
          if (isIsotopeShiftMatch(*it1, *it2, max_shift, c13_mass, mz_tol, mz_ppm))
          {
            // Delegate to Helper: Merge, annotate, log, and erase
            mergeAndAnnotateRescue(*it1, it2, candidate_map_index, out);
          }
        }
      }
    }

    postprocess_(maps, out);
  }

  void FeatureGroupingAlgorithmQT::group(const std::vector<FeatureMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmQT::group(const std::vector<ConsensusMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

} // namespace OpenMS