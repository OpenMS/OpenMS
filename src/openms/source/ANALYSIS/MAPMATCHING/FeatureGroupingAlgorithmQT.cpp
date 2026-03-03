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

#include <unordered_set>
#include <cmath>

using namespace std;

namespace OpenMS
{
  namespace
  {
    using ::OpenMS::ConsensusFeature;
    using ::OpenMS::Size;

    inline std::unordered_set<Size> buildMapIndexSet(const ConsensusFeature& cf)
    {
      std::unordered_set<Size> s;
      for (const auto& f : cf.getFeatures())
      {
        s.insert(f.getMapIndex());
      }
      return s;
    }

    // Compare an expected shifted m/z to two mz values, honoring ppm vs Da tolerance
    inline bool mzMatches(double mz1, double mz2, double expected_shift, bool mz_ppm, double mz_tol)
    {
      double diff = std::abs(mz1 - mz2 - expected_shift);
      if (mz_ppm)
      {
        // preventive guard: avoid division by zero (mz1 should not be zero for real peaks)
        if (mz1 <= 0.0) return false;
        double ppm = (diff / mz1) * 1e6;
        return ppm <= mz_tol;
      }
      else
      {
        return diff <= mz_tol;
      }
    }

    inline bool tryIsotopeShiftMatch(const ConsensusFeature& f1,
                                     const ConsensusFeature& f2,
                                     int max_shift,
                                     int max_assumed_charge,
                                     double c13_mass,
                                     bool mz_ppm,
                                     double mz_tol)
    {
      int charge1 = f1.getCharge();
      int charge2 = f2.getCharge();

      if (charge1 != 0 && charge2 != 0 && charge1 != charge2) return false;

      double mz1 = f1.getMZ();
      double mz2 = f2.getMZ();

      for (int k = 1; k <= max_shift; ++k)
      {
        if (charge1 == 0 && charge2 == 0)
        {
          // both unknown: try plausible charges
          for (int z = 1; z <= max_assumed_charge; ++z)
          {
            double expected_shift = (k * c13_mass) / static_cast<double>(z);
            if (mzMatches(mz1, mz2, expected_shift, mz_ppm, mz_tol)) return true;
          }
        }
        else
        {
          // at least one known: use the known one
          int z = (charge1 != 0) ? std::abs(charge1) : std::abs(charge2);
          double expected_shift = (k * c13_mass) / static_cast<double>(z);
          if (mzMatches(mz1, mz2, expected_shift, mz_ppm, mz_tol)) return true;
        }
      }

      return false;
    }
  } // unnamed namespace

  FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmQT");
    defaults_.insert("", QTClusterFinder().getParameters());

    defaults_.setValue("enable_isotope_shift_fallback", "false",
      "If true, performs a second linking pass to rescue features missing a monoisotopic peak.");
    defaults_.setValidStrings("enable_isotope_shift_fallback", {"true","false"});

    defaults_.setValue("max_isotope_shift", 2,
      "Maximum number of isotopic peaks to shift when matching missing M+0 features (e.g., 1 for M+1, 2 for M+2).");
    defaults_.setMinInt("max_isotope_shift", 1);

    defaults_.setValue("max_assumed_charge", 3,
      "Maximum charge state assumed when both features have unknown charge in isotope rescue.");
    defaults_.setMinInt("max_assumed_charge", 1);

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT() = default;

  template <typename MapType>
  void FeatureGroupingAlgorithmQT::group_(const vector<MapType>& maps, ConsensusMap& out)
  {
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
        OPENMS_PRETTY_FUNCTION, "At least two maps must be given!");
    }

    QTClusterFinder cluster_finder;
    cluster_finder.setParameters(param_.copy("", true));

    // PASS 1: Standard exact m/z clustering
    cluster_finder.run(maps, out);

    // PASS 2: Isotopic shift rescue (refactored into helpers to reduce method complexity)
    if (param_.getValue("enable_isotope_shift_fallback").toBool())
    {
      const Int max_shift = (Int)param_.getValue("max_isotope_shift");
      const Int max_assumed_charge = (Int)param_.getValue("max_assumed_charge");

      // Reuse QT tolerances (no hardcoding)
      const double rt_tol = param_.getValue("distance_RT:max_difference");
      const double mz_tol = param_.getValue("distance_MZ:max_difference");
      const bool mz_ppm = param_.getValue("distance_MZ:unit").toString() == "ppm";

      const double c13_mass = 1.0033548378;

      // Greedy pass to merge unlinked singleton consensus features
      for (auto it1 = out.begin(); it1 != out.end(); ++it1)
      {
        if (it1->getFeatures().size() != 1) continue;

        const auto maps_it1 = buildMapIndexSet(*it1);

        for (auto it2 = it1 + 1; it2 != out.end(); ++it2)
        {
          if (it2->getFeatures().size() != 1) continue;

          // Skip if any map index overlaps (robust multi-map guard)
          bool overlap = false;
          for (const auto& f : it2->getFeatures())
          {
            if (maps_it1.count(f.getMapIndex()))
            {
              overlap = true;
              break;
            }
          }
          if (overlap) continue;
          if (std::abs(it1->getRT() - it2->getRT()) > rt_tol) continue;
          if (tryIsotopeShiftMatch(*it1, *it2, max_shift, max_assumed_charge, c13_mass, mz_ppm, mz_tol))
          {
            it1->insert(*it2->getFeatures().begin());
            it1->setMetaValue("isotope_shift_rescued", "true");
            it2 = out.erase(it2);
            --it2;
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