// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
      @brief Base class for all feature grouping algorithms

      These algorithms group corresponding features in one map or across maps.
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithm :
    public DefaultParamHandler
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithm();

    /// Destructor
    ~FeatureGroupingAlgorithm() override;

    ///Applies the algorithm. The features in the input @p maps are grouped and the output is written to the consensus map @p out
    virtual void group(const std::vector<FeatureMap > & maps, ConsensusMap & out) = 0;

    ///Applies the algorithm. The consensus features in the input @p maps are grouped and the output is written to the consensus map @p out
    /// Algorithms not supporting ConsensusMap input should simply not override this method,
    /// as the base implementation will forward the data to the FeatureMap version of group()
    virtual void group(const std::vector<ConsensusMap> & maps, ConsensusMap & out);

    /// Transfers subelements (grouped features) from input consensus maps to the result consensus map
    void transferSubelements(const std::vector<ConsensusMap> & maps, ConsensusMap & out) const;


protected:

    /// after grouping by the subclasses, postprocess unassigned IDs, protein IDs and sort results in a
    /// consistent way
    template<class MapType>
    void postprocess_(const std::vector<MapType>& maps, ConsensusMap& out)
    {
      // add protein IDs and unassigned peptide IDs to the result map here,
      // to keep the same order as the input maps (useful for output later):
      auto& newIDs = out.getUnassignedPeptideIdentifications();
      Size map_idx = 0;

      for (typename std::vector<MapType>::const_iterator map_it = maps.begin();
           map_it != maps.end(); ++map_it)
      {
        // add protein identifications to result map:
        out.getProteinIdentifications().insert(
            out.getProteinIdentifications().end(),
            map_it->getProteinIdentifications().begin(),
            map_it->getProteinIdentifications().end());

        // assign the map_index to unassigned PepIDs as well.
        // for the assigned ones, this has to be done in the subclass.
        for (const PeptideIdentification& pepID : map_it->getUnassignedPeptideIdentifications())
        {
          auto newPepID = pepID;
          // Note: during linking of _consensus_Maps we have the problem that old identifications
          // should already have a map_index associated. Since we group the consensusFeatures only anyway
          // (without keeping the subfeatures) the method for now is to "re"-index based on the input file/map index.
          // Subfeatures have to be transferred in postprocessing if required
          // (see FeatureGroupingAlgorithm::transferSubelements as used in the TOPP tools, i.e. FeatureLinkerBase),
          // which also takes care of a re-re-indexing if the old map_index of the IDs was saved.
          newPepID.setMetaValue("map_index", map_idx);
          newIDs.push_back(newPepID);
        }
        map_idx++;
      }

      // canonical ordering for checking the results:
      out.sortByQuality();
      out.sortByMaps();
      out.sortBySize();
    }
private:
    ///Copy constructor is not implemented -> private
    FeatureGroupingAlgorithm(const FeatureGroupingAlgorithm &);
    ///Assignment operator is not implemented -> private
    FeatureGroupingAlgorithm & operator=(const FeatureGroupingAlgorithm &);



  };

} // namespace OpenMS

