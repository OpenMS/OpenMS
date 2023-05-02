// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
    for (const auto& cons : consensus)
    {
      DPosition<2> point(cons.getRT(), cons.getMZ());
      // for each sub-map in the consensus map, add a seed at the position of
      // this consensus feature:
      for (const auto& file : consensus.getColumnHeaders())
      {
        seed_lists[file.first].push_back(point);
      }
      // for each feature contained in the consensus feature, remove the seed of
      // the corresponding map:
      for (const auto& feat : cons.getFeatures())
      {
        seed_lists[feat.getMapIndex()].pop_back();
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
    for (const auto& seed : seeds)
    {
      Feature feature;
      feature.setRT(seed.getX());
      feature.setMZ(seed.getY());
      feature.setUniqueId(counter);
      features.push_back(feature);
      ++counter;
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
