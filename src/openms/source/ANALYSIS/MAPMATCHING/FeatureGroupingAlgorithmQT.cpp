// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmQT");
    defaults_.insert("", QTClusterFinder().getParameters());
    defaultsToParam_();
  }

  FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT()
  {
  }

  template <typename MapType>
  void FeatureGroupingAlgorithmQT::group_(const vector<MapType> & maps,
                                          ConsensusMap & out)
  {
    // check that the number of maps is ok:
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "At least two maps must be given!");
    }

    QTClusterFinder cluster_finder;
    cluster_finder.setParameters(param_.copy("", true));

    cluster_finder.run(maps, out);

    // add protein IDs and unassigned peptide IDs to the result map here,
    // to keep the same order as the input maps (useful for output later):
    for (typename vector<MapType>::const_iterator map_it = maps.begin();
         map_it != maps.end(); ++map_it)
    {
      // add protein identifications to result map:
      out.getProteinIdentifications().insert(
        out.getProteinIdentifications().end(),
        map_it->getProteinIdentifications().begin(),
        map_it->getProteinIdentifications().end());

      // add unassigned peptide identifications to result map:
      out.getUnassignedPeptideIdentifications().insert(
        out.getUnassignedPeptideIdentifications().end(),
        map_it->getUnassignedPeptideIdentifications().begin(),
        map_it->getUnassignedPeptideIdentifications().end());
    }

    // canonical ordering for checking the results:
    out.sortByQuality();
    out.sortByMaps();
    out.sortBySize();
    return;
  }

  void FeatureGroupingAlgorithmQT::group(const std::vector<FeatureMap> & maps,
                                         ConsensusMap & out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmQT::group(const std::vector<ConsensusMap> & maps,
                                         ConsensusMap & out)
  {
    group_(maps, out);
  }

} // namespace OpenMS
