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
#include <OpenMS/CONCEPT/Constants.h>

using namespace std;

namespace OpenMS
{

  namespace
  {
    // Helper 1: Charge compatibility check (Simplified per maintainer review)
    bool isChargeCompatible(const ConsensusFeature& a, const ConsensusFeature& b, Int& effective_charge)
    {
      const Int& charge_1 = a.getCharge();
      const Int& charge_2 = b.getCharge();

      // Handle both zero early
      if (charge_1 == 0 && charge_2 == 0) return false;

      // Charges must be identical if both are known (no std::abs needed)
      if (charge_1 != 0 && charge_2 != 0 && charge_1 != charge_2) return false;

      effective_charge = (charge_1 != 0) ? charge_1 : charge_2;
      return true;
    }

    // Helper 2: Isotope shift math and tolerance check using OpenMS constants
    bool isIsotopeShiftMatch(const ConsensusFeature& a, const ConsensusFeature& b, const std::vector<int>& allowed_shifts, Int effective_charge, double mz_tol, bool mz_ppm, int& matched_error)
    {
      double mz_diff = std::abs(a.getMZ() - b.getMZ());

      for (int k : allowed_shifts)
      {
        double expected_mz_shift = (k * Constants::C13C12_MASSDIFF_U) / effective_charge;
        double shift_diff = std::abs(mz_diff - expected_mz_shift);
        double tolerance = mz_ppm ? (mz_tol * a.getMZ() / 1e6) : mz_tol;

        if (shift_diff <= tolerance)
        {
          matched_error = k;
          return true;
        }
      }
      return false;
    }
  }

  FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmQT");
    defaults_.insert("", QTClusterFinder().getParameters());

    // Consolidated parameters into a single IntList
    defaults_.setValue("add_isotope_error", std::vector<int>(), "If not empty, performs a second linking pass to rescue features missing a monoisotopic peak. Specifies the allowed isotope errors (e.g., [1, 2] for M+1 and M+2).");

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT() = default;

  void FeatureGroupingAlgorithmQT::updateMembers_()
  {
    FeatureGroupingAlgorithm::updateMembers_();

    // Read params BEFORE running cluster finder
    add_isotope_error_ = param_.getValue("add_isotope_error");
    rt_tolerance_ = param_.getValue("distance_RT:max_difference");
    mz_tolerance_ = param_.getValue("distance_MZ:max_difference");
    mz_measure_is_ppm_ = param_.getValue("distance_MZ:unit").toString() == "ppm";
  }

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
    if (!add_isotope_error_.empty())
    {
      // 1. Sort by RT to allow for a high-performance O(N log N) sliding window
      out.sortByRT();

      for (auto it1 = out.begin(); it1 != out.end(); ++it1)
      {
        // Skip features that were already merged or aren't singletons
        if (it1->metaValueExists("merged_in_rescue") || it1->getFeatures().size() != 1) continue;

        for (auto it2 = it1 + 1; it2 != out.end(); ++it2)
        {
          if (it2->metaValueExists("merged_in_rescue") || it2->getFeatures().size() != 1) continue;

          // Sliding Window: Break inner loop early if RT difference exceeds tolerance
          if ((it2->getRT() - it1->getRT()) > rt_tolerance_)
          {
            break;
          }

          // Ensure candidate map is not already represented in the merged group
          const auto candidate_map_index = it2->getFeatures().begin()->getMapIndex();
          const bool map_already_present = std::any_of(
            it1->getFeatures().begin(), it1->getFeatures().end(),
            [candidate_map_index](const auto& handle)
            {
              return handle.getMapIndex() == candidate_map_index;
            });

          if (map_already_present) continue;

          Int effective_charge = 0;
          if (!isChargeCompatible(*it1, *it2, effective_charge)) continue;

          int matched_error = 0;
          if (isIsotopeShiftMatch(*it1, *it2, add_isotope_error_, effective_charge, mz_tolerance_, mz_measure_is_ppm_, matched_error))
          {
            // Merge the feature
            it1->insert(*it2->getFeatures().begin());

            // Annotate the integer isotope error per maintainer request
            it1->setMetaValue("isotope_error", matched_error);

            OPENMS_LOG_DEBUG << "FeatureGroupingAlgorithmQT: Isotope-shift rescue applied. Merged map "
                             << candidate_map_index << " into consensus feature.\n";

            // Mark it2 for safe O(N) deletion later
            it2->setMetaValue("merged_in_rescue", "true");
          }
        }
      }

      // Cleanup: Safely remove all singletons that were merged
      out.erase(std::remove_if(out.begin(), out.end(),
                [](const ConsensusFeature& f) { return f.metaValueExists("merged_in_rescue"); }),
                out.end());
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