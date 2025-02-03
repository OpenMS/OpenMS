// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/SeedListGenerator.h>
#include <map>

using namespace std;

namespace OpenMS
{
  SeedListGenerator::SeedListGenerator() = default;

  void SeedListGenerator::generateSeedList(const PeakMap& experiment,
                                           SeedList& seeds)
  {
    seeds.clear();
    for (PeakMap::ConstIterator exp_it = experiment.begin();
         exp_it != experiment.end(); ++exp_it)
    {
      if (exp_it->getMSLevel() == 2) // MS2 spectrum -> look for precursor
      {
        PeakMap::ConstIterator prec_it =
          experiment.getPrecursorSpectrum(exp_it);
        const vector<Precursor>& precursors = exp_it->getPrecursors();
        DPosition<2> point(prec_it->getRT(), precursors[0].getMZ());
        seeds.push_back(point);
      }
    }
  }

  void SeedListGenerator::generateSeedList(vector<PeptideIdentification>&
                                           peptides, SeedList& seeds,
                                           bool use_peptide_mass)
  {
    seeds.clear();
    for (PeptideIdentification& pep : peptides)
    {
      double mz;
      if (!pep.getHits().empty() && use_peptide_mass)
      {
        pep.sort();
        const PeptideHit& hit = pep.getHits().front();
        Int charge = hit.getCharge();
        mz = hit.getSequence().getMZ(charge);
      }
      else
      {
        mz = pep.getMZ();
      }
      DPosition<2> point(pep.getRT(), mz);
      seeds.push_back(point);
    }
  }

  void SeedListGenerator::generateSeedLists(const ConsensusMap& consensus,
                                            std::map<UInt64, SeedList>& seed_lists)
  {
    seed_lists.clear();
    // iterate over all consensus features...
    for (ConsensusMap::ConstIterator cons_it = consensus.begin();
         cons_it != consensus.end(); ++cons_it)
    {
      DPosition<2> point(cons_it->getRT(), cons_it->getMZ());
      // for each sub-map in the consensus map, add a seed at the position of
      // this consensus feature:
      for (ConsensusMap::ColumnHeaders::const_iterator file_it =
             consensus.getColumnHeaders().begin(); file_it !=
           consensus.getColumnHeaders().end(); ++file_it)
        seed_lists[file_it->first].push_back(point);
      // for each feature contained in the consensus feature, remove the seed of
      // the corresponding map:
      for (ConsensusFeature::HandleSetType::const_iterator feat_it =
             cons_it->getFeatures().begin(); feat_it !=
           cons_it->getFeatures().end(); ++feat_it)
      {
        seed_lists[feat_it->getMapIndex()].pop_back();
      }
      // this leaves seeds for maps where no feature was found near the
      // consensus position
    }
  }

  void SeedListGenerator::convertSeedList(const SeedList& seeds,
                                          FeatureMap& features)
  {
    features.clear(true); // "true" should really be a default value here...
    Size counter = 0;
    for (SeedList::const_iterator seed_it = seeds.begin();
         seed_it != seeds.end(); ++seed_it, ++counter)
    {
      Feature feature;
      feature.setRT(seed_it->getX());
      feature.setMZ(seed_it->getY());
      feature.setUniqueId(counter);
      features.push_back(feature);
    }
    // // assign unique ids:
    // features.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  }

  void SeedListGenerator::convertSeedList(const FeatureMap& features,
                                          SeedList& seeds)
  {
    seeds.clear();
    for (const Feature& feat : features)
    {
      DPosition<2> point(feat.getRT(), feat.getMZ());
      seeds.push_back(point);
    }
  }

} // namespace OpenMS
