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

FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT() :
  FeatureGroupingAlgorithm()
{
  setName("FeatureGroupingAlgorithmQT");
  defaults_.insert("", QTClusterFinder().getParameters());

  defaults_.setValue("enable_isotope_shift_fallback", "false",
    "If true, performs a second linking pass to rescue features missing a monoisotopic peak.");
  defaults_.setValidStrings("enable_isotope_shift_fallback", {"true","false"});

  defaults_.setValue("max_isotope_shift", 2,
    "Maximum number of isotopic peaks to shift when matching missing M+0 features.");
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

  // PASS 1: Standard QT clustering
  cluster_finder.run(maps, out);

  // PASS 2: Optional isotope rescue
  if (param_.getValue("enable_isotope_shift_fallback").toBool())
  {
    rescueIsotopeShift_(out);
  }

  postprocess_(maps, out);
}

void FeatureGroupingAlgorithmQT::rescueIsotopeShift_(ConsensusMap& out)
{
  const Int max_shift = (Int)param_.getValue("max_isotope_shift");
  const Int max_assumed_charge = (Int)param_.getValue("max_assumed_charge");

  const double rt_tol = param_.getValue("distance_RT:max_difference");
  const double mz_tol = param_.getValue("distance_MZ:max_difference");
  const bool mz_ppm = param_.getValue("distance_MZ:unit").toString() == "ppm";

  const double c13_mass = 1.0033548378;

  auto mzMatches = [&](double mz1, double mz2, double expected_shift)
  {
    double diff = std::abs(mz1 - mz2 - expected_shift);

    if (mz_ppm)
    {
      if (mz1 <= 0.0) return false;
      double ppm = (diff / mz1) * 1e6;
      return ppm <= mz_tol;
    }
    return diff <= mz_tol;
  };

  for (auto it1 = out.begin(); it1 != out.end(); ++it1)
  {
    if (it1->getFeatures().size() != 1) continue;

    // Collect map indices of first group
    std::unordered_set<Size> maps_it1;
    for (const auto& f : it1->getFeatures())
    {
      maps_it1.insert(f.getMapIndex());
    }

    for (auto it2 = it1 + 1; it2 != out.end(); ++it2)
    {
      if (it2->getFeatures().size() != 1) continue;

      // Robust map uniqueness guard
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

      // RT tolerance check
      if (std::abs(it1->getRT() - it2->getRT()) > rt_tol) continue;

      int charge1 = it1->getCharge();
      int charge2 = it2->getCharge();

      if (charge1 != 0 && charge2 != 0 && charge1 != charge2)
        continue;

      bool match_found = false;

      for (int k = 1; k <= max_shift && !match_found; ++k)
      {
        if (charge1 == 0 && charge2 == 0)
        {
          // Try plausible charges
          for (int z = 1; z <= max_assumed_charge && !match_found; ++z)
          {
            double expected_shift = (k * c13_mass) / (double)z;
            if (mzMatches(it1->getMZ(), it2->getMZ(), expected_shift))
              match_found = true;
          }
        }
        else
        {
          int z = (charge1 != 0) ? std::abs(charge1) : std::abs(charge2);
          double expected_shift = (k * c13_mass) / (double)z;

          if (mzMatches(it1->getMZ(), it2->getMZ(), expected_shift))
            match_found = true;
        }
      }

      if (match_found)
      {
        it1->insert(*it2->getFeatures().begin());
        it1->setMetaValue("isotope_shift_rescued", "true");

        it2 = out.erase(it2);
        --it2;
      }
    }
  }
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