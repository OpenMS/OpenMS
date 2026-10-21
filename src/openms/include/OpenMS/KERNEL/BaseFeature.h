// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <optional>

namespace OpenMS
{
  class FeatureHandle;

  /**
    @brief A basic LC-MS feature.

    This class represents a "minimal" feature, defined by a position in RT and m/z, intensity,
    charge, quality, and annotated peptides. Most classes dealing with features will use the
    subclasses Feature or ConsensusFeature directly. However, algorithms that rely on very general
    characteristics of features can use this class to provide a unified solution for both "normal"
    features and consensus features.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI BaseFeature : public RichPeak2D
  {
public:
    /// @name Type definitions
    ///@{
    /// Type of quality values
    typedef float QualityType;
    /// Type of charge values
    typedef Int ChargeType;
    /// Type of feature width/FWHM (RT)
    typedef float WidthType;

    /// state of identification, use getAnnotationState() to query it
    enum class AnnotationState
    {
      FEATURE_ID_NONE,
      FEATURE_ID_SINGLE,
      FEATURE_ID_MULTIPLE_SAME,
      FEATURE_ID_MULTIPLE_DIVERGENT,
      SIZE_OF_ANNOTATIONSTATE
    };

    static const std::string NamesOfAnnotationState[static_cast<size_t>(AnnotationState::SIZE_OF_ANNOTATIONSTATE)];
    ///@}

    /// @name Constructors and Destructor
    ///@{
    /// Default constructor
    BaseFeature();

    /// Copy constructor
    BaseFeature(const BaseFeature& feature) = default;

    /// Move constructor
    /// Note: can't be "noexcept = default" because of missing noexcept on some standard containers
    /// so we need to explicitly define it noexcept and provide an implementation.
    BaseFeature(BaseFeature&& feature) noexcept
      : RichPeak2D(std::move(feature))
    {
      quality_ = feature.quality_;
      charge_ = feature.charge_;
      width_ = feature.width_;
      // Note: will terminate program if move assignment throws because of noexcept
      // but we can't recover in that case anyways and we need to mark it noexcept for the move.
      primary_id_ = std::move(feature.primary_id_);
      id_matches_ = std::move(feature.id_matches_);
      id_queries_ = std::move(feature.id_queries_);
    }

    /// Constructor from raw data point
    explicit BaseFeature(const Peak2D& point);

    /// Constructor from raw data point with meta information
    explicit BaseFeature(const RichPeak2D& point);

    /// Constructor from a featurehandle
    explicit BaseFeature(const FeatureHandle& fh);

    /// Destructor
    ~BaseFeature() override;
    ///@}

    /// @name Quality methods
    ///@{
    /// Non-mutable access to the overall quality
    QualityType getQuality() const;
    /// Set the overall quality
    void setQuality(QualityType q);
    /// Compare by quality
    struct QualityLess
    {
      bool operator()(const BaseFeature& left, const BaseFeature& right) const
      {
        return left.getQuality() < right.getQuality();
      }

      bool operator()(const BaseFeature& left, const QualityType& right) const
      {
        return left.getQuality() < right;
      }

      bool operator()(const QualityType& left, const BaseFeature& right) const
      {
        return left < right.getQuality();
      }

      bool operator()(const QualityType& left, const QualityType& right) const
      {
        return left < right;
      }

    };
    ///@}

    /// Non-mutable access to the features width (full width at half max, FWHM)
    WidthType getWidth() const;
    /// Set the width of the feature (FWHM)
    void setWidth(WidthType fwhm);

    /// Non-mutable access to charge state
    const ChargeType& getCharge() const;

    /// Set charge state
    void setCharge(const ChargeType& ch);

    /// Assignment operator
    BaseFeature& operator=(const BaseFeature& rhs) = default;

    /// Move Assignment operator
    BaseFeature& operator=(BaseFeature&& rhs) & = default;

    /// Equality operator
    bool operator==(const BaseFeature& rhs) const;

    /// Inequality operator
    bool operator!=(const BaseFeature& rhs) const;

    /**
      @brief State of the identifications of @p peptides (e.g. the peptide identifications of a feature in a format, see
      MapPeptideIdentifications): none, one, several of the same sequence or of different sequences (by the best hit of
      each identification)
    */
    static AnnotationState annotationState(const PeptideIdentificationList& peptides);
    /**
      @brief State of the identifications this feature links in @p data, as annotationState() gives it for their
      peptide identifications: none, one, several of the same molecule or of different molecules (by the best linked
      match of each identification)

      @throws Exception::MissingInformation if a link refers to an identification or match that @p data does not have
    */
    AnnotationState getAnnotationState(const IdentificationData& data) const;

    /// @name Functions for dealing with identifications in new format
    ///@{
    /// has a primary ID (peptide, RNA, compound) been assigned?
    bool hasPrimaryID() const;

    /**
       @brief Return the primary ID (peptide, RNA, compound) assigned to this feature.

       @throw Exception::MissingInformation if no ID was assigned
    */
    const IdentificationData::MoleculeIdentity& getPrimaryID() const;

    /// clear any primary ID that was assigned
    void clearPrimaryID();

    /// set the primary ID (peptide, RNA, compound) for this feature
    void setPrimaryID(const IdentificationData::MoleculeIdentity& id);

    /// immutable access to the set of matches (e.g. PSMs) with IDs for this feature
    const std::set<IdentificationData::MatchReference>& getIDMatches() const;

    /// mutable access to the set of matches (e.g. PSMs) with IDs for this feature
    std::set<IdentificationData::MatchReference>& getIDMatches();

    const std::set<IdentificationData::QueryReference>& getIDQueries() const
    { return id_queries_; }
    std::set<IdentificationData::QueryReference>& getIDQueries()
    { return id_queries_; }
    void addIDQuery(const IdentificationData::QueryReference& reference)
    { id_queries_.insert(reference); }

    /**
       @brief The identifications this feature links, directly or through their matches

       @throw Exception::MissingInformation if a match link refers to a match that @p data does not have
    */
    std::set<IdentificationData::QueryReference> getLinkedIDQueries(const IdentificationData& data) const;

    /**
       @brief The identifications this feature links, each with its linked matches (its peptide identifications in a format)

       In the order of the peptide identifications that an export of @p data gives the feature.

       @throw Exception::MissingInformation if a link refers to something that @p data does not have
    */
    std::vector<IdentificationData::QueryMatches> getLinkedIdentifications(const IdentificationData& data) const;

    /**
       @brief The linked match with the best primary score, with its identification (the best hit of the peptide
       identifications of the feature)

       Of matches that score equally, the first in the order of getLinkedIdentifications() is taken.
       Matches without a primary score value are skipped.

       @return An entry with exactly one match, or nothing if the feature links no match with a score
       @throw Exception::MissingInformation if a link refers to something that @p data does not have
    */
    std::optional<IdentificationData::QueryMatches> getBestLinkedMatch(const IdentificationData& data) const;

    /// add an ID match (e.g. PSM) for this feature
    void addIDMatch(IdentificationData::MatchReference ref);

    ///@}

protected:

    /// Overall quality measure of the feature
    QualityType quality_;

    /// Charge of the peptide represented by this feature.  The default value is 0, which represents an unknown charge state.
    ChargeType charge_;

    /// Width (FWHM) for the feature. The default value is 0.0, a feature finding algorithm can compute this form the model.
    WidthType width_;

    /// primary ID (peptide, RNA, compound) assigned to this feature
    std::optional<IdentificationData::MoleculeIdentity> primary_id_;

    /// set of observation matches (e.g. PSMs) with IDs for this feature
    std::set<IdentificationData::MatchReference> id_matches_;
    std::set<IdentificationData::QueryReference> id_queries_;
  };

} // namespace OpenMS
