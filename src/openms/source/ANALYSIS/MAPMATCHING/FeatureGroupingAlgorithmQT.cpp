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

using namespace std;

namespace OpenMS
{

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
      Int max_shift = (Int)param_.getValue("max_isotope_shift");
      double c13_mass = 1.0033548; // Mass difference between 13C and 12C

      // Greedy pass to merge unlinked features (singletons)
      for (auto it1 = out.begin(); it1 != out.end(); ++it1)
      {
          // Only look at features that failed to link in Pass 1
          if (it1->getFeatures().size() != 1) continue;

          for (auto it2 = it1 + 1; it2 != out.end(); ++it2)
          {
              if (it2->getFeatures().size() != 1) continue;

              // Ensure they come from different maps
              if (it1->getFeatures().begin()->getMapIndex() == it2->getFeatures().begin()->getMapIndex()) continue;

              // Basic RT check (assuming a generous 20s tolerance for fallback)
              if (std::abs(it1->getRT() - it2->getRT()) > 20.0) continue;

              // Charges must match (if known)
              int charge = it1->getCharge();
              if (charge == 0) charge = 1; // Prevent division by zero
              if (charge != it2->getCharge() && it2->getCharge() != 0) continue;

              double mz_diff = std::abs(it1->getMZ() - it2->getMZ());
              bool match_found = false;

              // Check for isotope shifts (M+1, M+2, etc.)
              for (int k = 1; k <= max_shift; ++k)
              {
                  double expected_shift = (k * c13_mass) / std::abs(charge);

                  // 0.05 Da tolerance window for the shifted m/z
                  if (std::abs(mz_diff - expected_shift) < 0.05)
                  {
                      match_found = true;
                      break;
                  }
              }

              if (match_found)
              {
                  // 1. Merge the feature from it2 into it1
                  it1->insert(*it2->getFeatures().begin());

                  // 2. Annotate the rescue for downstream data integrity
                  it1->setMetaValue("isotope_shift_rescued", "true");

                  // 3. Remove it2 from the map so it isn't processed again
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
