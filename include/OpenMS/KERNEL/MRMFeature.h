// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MRMFEATURE_H
#define OPENMS_KERNEL_MRMFEATURE_H

#include <OpenMS/KERNEL/Feature.h>

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
    /// Peak group score type
    typedef std::map<String, double> PGScoresType;
    //@}

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    MRMFeature();

    /// Copy constructor
    MRMFeature(const MRMFeature &rhs);

    /// Assignment operator
    MRMFeature & operator=(const MRMFeature & rhs);

    /// Destructor
    virtual ~MRMFeature();
    //@}

    ///@name Accessors
    //@{
    /// get all peakgroup scores
    const PGScoresType & getScores() const;

    /// get a single peakgroup score
    double getScore(const String & score_name);

    /// get a specified feature
    Feature & getFeature(String key);

    /// set all peakgroup scores
    void setScores(const PGScoresType & scores);

    /// set a single peakgroup score
    void addScore(const String & score_name, double score);

    /// Adds an feature from a single chromatogram into the feature.
    void addFeature(Feature & feature, String key);

    const std::vector<Feature> & getFeatures() const;

    void getFeatureIDs(std::vector<String> & result) const;
    //@}

protected:

    FeatureListType features_;

    /// peak group scores
    PGScoresType pg_scores_;

    /// map native ids to the features
    std::map<String, int> feature_map_;

  };
}

#endif
