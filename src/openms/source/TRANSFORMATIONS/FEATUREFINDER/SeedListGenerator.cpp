// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>

using namespace std;

namespace OpenMS
{
  SeedListGenerator::SeedListGenerator()
  {
  }

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
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      double mz;
      if (!pep_it->getHits().empty() && use_peptide_mass)
      {
        pep_it->sort();
        const PeptideHit& hit = pep_it->getHits().front();
        Int charge = hit.getCharge();
        mz = hit.getSequence().getMonoWeight(Residue::Full, charge) /
             double(charge);
      }
      else
      {
        mz = pep_it->getMZ();
      }
      DPosition<2> point(pep_it->getRT(), mz);
      seeds.push_back(point);
    }
  }

  void SeedListGenerator::generateSeedLists(const ConsensusMap& consensus,
                                            Map<UInt64, SeedList>& seed_lists)
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
    for (FeatureMap::ConstIterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      DPosition<2> point(feat_it->getRT(), feat_it->getMZ());
      seeds.push_back(point);
    }
  }

} // namespace OpenMS
