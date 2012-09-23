// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>


namespace OpenMS
{

  FeatureGroupingAlgorithmUnlabeled::FeatureGroupingAlgorithmUnlabeled() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmUnlabeled");
    defaults_.insert("", StablePairFinder().getParameters());
    defaultsToParam_();
  }

  FeatureGroupingAlgorithmUnlabeled::~FeatureGroupingAlgorithmUnlabeled()
  {
  }

  void FeatureGroupingAlgorithmUnlabeled::group(const std::vector<FeatureMap<> > & maps, ConsensusMap & out)
  {
    // check that the number of maps is ok
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "At least two maps must be given!");
    }

    // define reference map (the one with most peaks)
    Size reference_map_index = 0;
    Size max_count = 0;
    for (Size m = 0; m < maps.size(); ++m)
    {
      if (maps[m].size() > max_count)
      {
        max_count = maps[m].size();
        reference_map_index = m;
      }
    }

    std::vector<ConsensusMap> input(2);

    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    ConsensusMap::convert(reference_map_index, maps[reference_map_index],
                          input[0]);

    // loop over all other maps, extend the groups
    StablePairFinder pair_finder;
    pair_finder.setParameters(param_.copy("", true));

    for (Size i = 0; i < maps.size(); ++i)
    {
      if (i != reference_map_index)
      {
        ConsensusMap::convert(i, maps[i], input[1]);
        // compute the consensus of the reference map and map i
        ConsensusMap result;
        pair_finder.run(input, result);
        input[0].swap(result);
      }
    }

    // replace result with temporary map
    out.swap(input[0]);
    // copy back the input maps (they have been deleted while swapping)
    out.getFileDescriptions() = input[0].getFileDescriptions();

    // add protein IDs and unassigned peptide IDs to the result map here,
    // to keep the same order as the input maps (useful for output later)
    for (std::vector<FeatureMap<> >::const_iterator map_it = maps.begin();
         map_it != maps.end(); ++map_it)
    {
      // add protein identifications to result map
      out.getProteinIdentifications().insert(
        out.getProteinIdentifications().end(),
        map_it->getProteinIdentifications().begin(),
        map_it->getProteinIdentifications().end());

      // add unassigned peptide identifications to result map
      out.getUnassignedPeptideIdentifications().insert(
        out.getUnassignedPeptideIdentifications().end(),
        map_it->getUnassignedPeptideIdentifications().begin(),
        map_it->getUnassignedPeptideIdentifications().end());
    }

    // canonical ordering for checking the results, and the ids have no real meaning anyway
#if 1 // the way this was done in DelaunayPairFinder and StablePairFinder
    out.sortByMZ();
#else
    out.sortByQuality();
    out.sortByMaps();
    out.sortBySize();
#endif

    return;
  }

} // namespace OpenMS
