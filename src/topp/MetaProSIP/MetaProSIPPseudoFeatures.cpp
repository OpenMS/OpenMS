// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MetaProSIPPseudoFeatures.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/Peak2D.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

namespace OpenMS
{
Size MetaProSIPPseudoFeatures::addUnassignedIdentifications(FeatureMap& features, const std::string& feature_type)
{
  std::vector<Feature> pseudo_features;
  for (const auto& unassigned : features.getUnassignedIdentifications())
  {
    if (unassigned.matches.empty()) { continue; }
    Feature f;
    f.setMetaValue("feature_type", feature_type);
    f.setRT(unassigned.query->rt.value_or(0.0));
    // take sequence of first match to calculate ground truth mz
    const IdentificationData::Match& first = *unassigned.matches[0];
    Int charge = first.charge;
    if (charge == 0) { continue; }
    f.setMZ(AASequence::fromString(first.representation).getMZ(charge));
    // the pseudo feature links the identification (with its matches)
    for (const auto* match : unassigned.matches)
    {
      f.addIDMatch({unassigned.run->getUuid(), match->getId()});
    }
    pseudo_features.push_back(f);
  }
  for (auto& f : pseudo_features)
  {
    features.push_back(std::move(f));
  }
  features.updateRanges();
  return pseudo_features.size();
}

Size MetaProSIPPseudoFeatures::addUnidentifiedSpectra(FeatureMap& features, const PeakMap& ms2, const std::string& feature_type)
{
  // extract rt and mz of all identified precursors and store them in blacklist
  std::vector<Peak2D> blacklisted_precursors;
  const auto blacklist = [&blacklisted_precursors](const IdentificationData::QueryMatches& identification)
  {
    if (identification.matches.empty())
    {
      return;
    }
    // Feature with id found so we don't need to generate averagine id. Find MS2 in experiment and blacklist it.
    Peak2D p;
    p.setRT(identification.query->rt.value_or(0.0));
    p.setMZ(identification.query->mz.value_or(0.0));
    blacklisted_precursors.push_back(p);
  };
  // in features
  for (const Feature& feature : features) // for each peptide feature
  {
    for (const auto& identification : feature.getLinkedIdentifications(features.getIdentificationData()))
    {
      blacklist(identification);
    }
  }
  // and in unassigned ids
  for (const auto& identification : features.getUnassignedIdentifications())
  {
    blacklist(identification);
  }

  // find index of all precursors that have been blacklisted
  std::vector<Size> blacklist_idx;
  for (const Peak2D& precursor : blacklisted_precursors)
  {
    PeakMap::const_iterator map_rt_begin = ms2.RTBegin(-std::numeric_limits<double>::max());
    PeakMap::const_iterator rt_begin = ms2.RTBegin(precursor.getRT() - 1e-5);
    blacklist_idx.push_back(std::distance(map_rt_begin, rt_begin));
  }

  Size added = 0;
  for (Size i = 0; i != ms2.size(); ++i)
  {
    // precursor not blacklisted?
    if (std::find(blacklist_idx.begin(), blacklist_idx.end(), i) == blacklist_idx.end() && !ms2[i].getPrecursors().empty())
    {
      // a pseudo feature for the unidentified spectrum, modelled as an averagine peptide
      Feature f;
      // the averagine model takes the place of an identified peptide: no sequence, the precursor's charge
      f.setCharge(ms2[i].getPrecursors()[0].getCharge());
      f.setRT(ms2[i].getRT());
      f.setMZ(ms2[i].getPrecursors()[0].getMZ());
      f.setMetaValue("feature_type", feature_type);
      features.push_back(f);
      ++added;
    }
  }
  features.updateRanges();
  return added;
}
} // namespace OpenMS
