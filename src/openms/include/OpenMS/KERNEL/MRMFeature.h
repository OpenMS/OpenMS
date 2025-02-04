// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>

namespace OpenMS
{
  /**
    @brief A multi-chromatogram MRM feature 

    An MRM feature contains corresponding features in multiple chromatograms,
    it is thus a representation of a peak group. The individual features in
    each chromatogram are represented by OpenMS Features.

  */
  class OPENMS_DLLAPI MRMFeature :
    public Feature
  {
public:

    ///Type definitions
    //@{
    /// Feature list type
    typedef std::vector<Feature> FeatureListType;
    //@}

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    MRMFeature();

    /// Copy constructor
    MRMFeature(const MRMFeature &rhs);

    /// Move constructor
    MRMFeature(MRMFeature &&rhs) = default;

    /// Assignment operator
    MRMFeature & operator=(const MRMFeature & rhs);

    /// Move assignment operator
    MRMFeature& operator=(MRMFeature&&) & = default;

    /// Destructor
    ~MRMFeature() override;
    //@}

    ///@name Accessors
    //@{

    /// get all peakgroup scores
    const OpenSwath_Scores & getScores() const;

    /// get all peakgroup scores
    OpenSwath_Scores & getScores();

    /// get a specified feature
    Feature & getFeature(const String& key);

    /// get a specified feature (const)
    const Feature & getFeature(const String& key) const;

    /// set all peakgroup scores
    void setScores(const OpenSwath_Scores & scores);

    /// set a single peakgroup score
    void addScore(const String & score_name, double score);

    /// Adds an feature from a single chromatogram into the feature.
    void addFeature(const Feature & feature, const String& key);

    void addFeature(Feature && feature, const String& key);

    /// get a list of features
    const std::vector<Feature> & getFeatures() const;

    /// get a list of IDs of available features
    void getFeatureIDs(std::vector<String> & result) const;

    /// Adds a precursor feature from a single chromatogram into the feature.
    void addPrecursorFeature(const Feature & feature, const String& key);

    void addPrecursorFeature(Feature && feature, const String& key);

    /// get a list of IDs of available precursor features
    void getPrecursorFeatureIDs(std::vector<String> & result) const;

    /// get a specified precursor feature
    Feature & getPrecursorFeature(const String& key);

    /// get a specified precursor feature (const)
    const Feature & getPrecursorFeature(const String& key) const;

    void IDScoresAsMetaValue(bool decoy, const OpenSwath_Ind_Scores& idscores);
    //@}

protected:

    FeatureListType features_;

    FeatureListType precursor_features_;

    /// peak group scores
    OpenSwath_Scores pg_scores_;

    /// map native ids to the features
    std::map<String, int> feature_map_;

    /// map native ids to the precursor features
    std::map<String, int> precursor_feature_map_;

  };
}


