// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

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
    enum AnnotationState
    {
      FEATURE_ID_NONE,
      FEATURE_ID_SINGLE,
      FEATURE_ID_MULTIPLE_SAME,
      FEATURE_ID_MULTIPLE_DIVERGENT,
      SIZE_OF_ANNOTATIONSTATE
    };

    static const std::string NamesOfAnnotationState[SIZE_OF_ANNOTATIONSTATE];
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
      peptides_ = std::move(feature.peptides_);
      primary_id_ = std::move(feature.primary_id_);
      id_matches_ = std::move(feature.id_matches_);
    }

    /// Copy constructor with a new map_index
    BaseFeature(const BaseFeature& rhs, UInt64 map_index);

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

    /// @name Functions for dealing with identifications in legacy format
    ///@{
    /// returns a const reference to the PeptideIdentification vector
    const std::vector<PeptideIdentification>& getPeptideIdentifications() const;

    /// returns a mutable reference to the PeptideIdentification vector
    std::vector<PeptideIdentification>& getPeptideIdentifications();

    /// sets the PeptideIdentification vector
    void setPeptideIdentifications(const std::vector<PeptideIdentification>& peptides);

    /// sorts PeptideIdentifications, assuming they have the same scoreType.
    void sortPeptideIdentifications();
    ///@}

    /// state of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only
    AnnotationState getAnnotationState() const;

    /// @name Functions for dealing with identifications in new format
    ///@{
    /// has a primary ID (peptide, RNA, compound) been assigned?
    bool hasPrimaryID() const;

    /**
       @brief Return the primary ID (peptide, RNA, compound) assigned to this feature.

       @throw Exception::MissingInformation if no ID was assigned
    */
    const IdentificationData::IdentifiedMolecule& getPrimaryID() const;

    /// clear any primary ID that was assigned
    void clearPrimaryID();

    /// set the primary ID (peptide, RNA, compound) for this feature
    void setPrimaryID(const IdentificationData::IdentifiedMolecule& id);

    /// immutable access to the set of matches (e.g. PSMs) with IDs for this feature
    const std::set<IdentificationData::ObservationMatchRef>& getIDMatches() const;

    /// mutable access to the set of matches (e.g. PSMs) with IDs for this feature
    std::set<IdentificationData::ObservationMatchRef>& getIDMatches();

    /// add an ID match (e.g. PSM) for this feature
    void addIDMatch(IdentificationData::ObservationMatchRef ref);

    /*!
      @brief Update ID references (primary ID, matches) for this feature

      This is needed e.g. after the IdentificationData instance containing the referenced data has been copied.
    */
    void updateIDReferences(const IdentificationData::RefTranslator& trans);
    ///@}

protected:

    /// Overall quality measure of the feature
    QualityType quality_;

    /// Charge of the peptide represented by this feature.  The default value is 0, which represents an unknown charge state.
    ChargeType charge_;

    /// Width (FWHM) for the feature. The default value is 0.0, a feature finding algorithm can compute this form the model.
    WidthType width_;

    /// PeptideIdentifications belonging to the feature
    std::vector<PeptideIdentification> peptides_;

    /// primary ID (peptide, RNA, compound) assigned to this feature
    std::optional<IdentificationData::IdentifiedMolecule> primary_id_;

    /// set of observation matches (e.g. PSMs) with IDs for this feature
    std::set<IdentificationData::ObservationMatchRef> id_matches_;
  };

} // namespace OpenMS
