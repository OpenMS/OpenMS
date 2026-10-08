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

#include <functional>
#include <set>

namespace OpenMS
{

  /**
      @brief Base class for all feature grouping algorithms

      These algorithms group corresponding features in one map or across maps.

      The result has the identifications of the maps (see groupIdentifications()). Maps with peptide identifications
      are grouped on identification data converted from them, and the result gets peptide identifications back.
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

    /**
        @brief Transfers subelements (grouped features) from input consensus maps to the result consensus map

        The map indices of the identifications follow: an identification that had a map index in its input map (saved
        as meta value "old_map_index" before grouping) gets the index of that map in the result, the map index of
        others is removed.
    */
    void transferSubelements(const std::vector<ConsensusMap> & maps, ConsensusMap & out) const;

    /// The identifications of a map whose features are grouped, with the links of its features and their subordinates
    /// (see groupIdentifications())
    struct OPENMS_DLLAPI MapIdentifications
    {
      IdentificationData data;
      std::set<IdentificationData::QueryReference> queries;
      std::set<IdentificationData::MatchReference> matches;
    };

    /// The identifications of @p map, to give them to the map that groups its features (see groupIdentifications())
    static MapIdentifications getMapIdentifications(const FeatureMap& map);
    static MapIdentifications getMapIdentifications(const ConsensusMap& map);

    /**
        @brief Give @p grouped the identifications of the maps whose features it groups, the counterpart of the peptide
        identifications that grouping copies from the grouped features and adds as unassigned ones

        The identifications of the i-th map are marked with its index (meta value "map_index"), the map index of its
        features in @p grouped. Identifications and matches that features of the map link, but no feature of @p grouped
        (e.g. those of subordinates or of features left out), are dropped; unassigned ones are kept.

        The runs of the maps must be distinct (see IdentificationDataConverter::withIdentificationData()).
    */
    static void groupIdentifications(std::vector<MapIdentifications> maps, ConsensusMap& grouped);
    static void groupIdentifications(const std::vector<FeatureMap>& maps, ConsensusMap& grouped);
    static void groupIdentifications(const std::vector<ConsensusMap>& maps, ConsensusMap& grouped);

protected:

    /**
        @brief After grouping by the subclasses, give the result the identifications of the maps and sort it in a
        consistent way

        Maps with identification data give theirs (see groupIdentifications()); otherwise the protein identifications
        and the unassigned peptide identifications (with the map index) of the maps are added to the result.
    */
    void postprocess_(const std::vector<FeatureMap>& maps, ConsensusMap& out) const;
    void postprocess_(const std::vector<ConsensusMap>& maps, ConsensusMap& out) const;

    /**
        @brief Group the @p maps with their identifications as identification data

        Calls @p group with the maps prepared by IdentificationDataConverter::withIdentificationData(). If a map has
        peptide identifications, so does @p out afterwards (IdentificationDataConverter::exportConsensusIDs()).
    */
    static void groupWithIdentificationData_(const std::vector<FeatureMap>& maps, ConsensusMap& out,
                                             const std::function<void(const std::vector<FeatureMap>&)>& group);
    static void groupWithIdentificationData_(const std::vector<ConsensusMap>& maps, ConsensusMap& out,
                                             const std::function<void(const std::vector<ConsensusMap>&)>& group);
private:
    ///Copy constructor is not implemented -> private
    FeatureGroupingAlgorithm(const FeatureGroupingAlgorithm &);
    ///Assignment operator is not implemented -> private
    FeatureGroupingAlgorithm & operator=(const FeatureGroupingAlgorithm &);



  };

} // namespace OpenMS

