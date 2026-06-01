// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/KERNEL/BaseFeature.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <vector>

namespace OpenMS
{

  /** @brief An LC-MS feature.

  The Feature class is used to describe the two-dimensional signal caused by an
  analyte. It can store a charge state and a list of peptide identifications
  (for peptides). The area occupied by the Feature in the LC-MS data set is
  represented by a list of convex hulls (one for each isotopic peak). There is
  also a convex hull for the entire Feature. The model description can store
  the parameters of a two-dimensional theoretical model of the underlying
  signal in LC-MS. Currently, non-peptide compounds are also represented as
  features.

  By convention in %OpenMS, the position of a feature is defined as maximum
  position of the model for the retention time dimension and the mass of the
  monoisotopic peak for the m/z dimension. The intensity of a feature is
  (proportional to) its total ion count.

  Feature is derived from RichPeak2D. Also inherited is a MetaInfoInterface.
  Features as usually are contained in a FeatureMap. See also FeatureHandle and
  ConsensusFeature.

  @ingroup Kernel
  */
  class OPENMS_DLLAPI Feature :
    public BaseFeature
  {
public:
    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    Feature();

    /// explicit C'tor from BaseFeature
    explicit Feature(const BaseFeature& base);

    /// Copy constructor
    Feature(const Feature& feature);

    /// Move constructor
    Feature(Feature&&) noexcept;

    /// Destructor
    ~Feature() override;
    //@}

    /// @name Model and quality methods
    //@{
    /// Non-mutable access to the overall quality
    QualityType getOverallQuality() const;

    /// Set the overall quality
    void setOverallQuality(QualityType q);

    /// Non-mutable access to the quality in dimension c
    QualityType getQuality(Size index) const;
    /// Set the quality in dimension c
    void setQuality(Size index, QualityType q);

    /// Compare by quality
    typedef QualityLess OverallQualityLess;

    //@}

    ///@name Convex hulls and bounding box
    //@{
    /// Non-mutable access to the convex hulls
    const std::vector<ConvexHull2D>& getConvexHulls() const;
    /// Mutable access to the convex hulls of single mass traces
    std::vector<ConvexHull2D>& getConvexHulls();
    /// Set the convex hulls of single mass traces
    void setConvexHulls(const std::vector<ConvexHull2D>& hulls);

    /**
      @brief Returns the overall convex hull of the feature (calculated from the convex hulls of the mass traces)

      @note the bounding box of the feature can be accessed through the returned convex hull
    */
    ConvexHull2D& getConvexHull() const;

    /// Returns if the mass trace convex hulls of the feature enclose the position specified by @p rt and @p mz
    bool encloses(double rt, double mz) const;
    //@}

    /// Assignment operator
    Feature& operator=(const Feature& rhs);

    /// Move assignment operator
    Feature& operator=(Feature&&) & noexcept;

    /// Equality operator
    bool operator==(const Feature& rhs) const;

    /// immutable access to subordinate features
    const std::vector<Feature>& getSubordinates() const;

    /// mutable access to subordinate features
    std::vector<Feature>& getSubordinates();

    /// mutable access to subordinate features
    void setSubordinates(const std::vector<Feature>& rhs);

    /**
      @brief Applies a member function of Type to the feature (including subordinates).
      The returned values are accumulated.

      <b>Example:</b>  The following will print the number of features (parent feature and subordinates) with invalid unique ids:
      @code
      Feature f;
      (...)
      std::cout << f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
      @endcode
      See e.g. UniqueIdInterface for what else can be done this way.
    */
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)())
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (std::vector<Feature>::iterator iter = subordinates_.begin(); iter != subordinates_.end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    /// The "const" variant.
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)() const) const
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (std::vector<Feature>::const_iterator iter = subordinates_.begin(); iter != subordinates_.end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    /*!
      @brief Update ID references (primary ID, input matches) for this feature and any subfeatures

      This is needed e.g. after the IdentificationData instance containing the referenced data has been copied.
    */
    void updateAllIDReferences(const IdentificationData::RefTranslator& trans);

protected:

    /// Quality measures for each dimension
    QualityType qualities_[2];

    /// Array of convex hulls (one for each mass trace)
    std::vector<ConvexHull2D> convex_hulls_;

    /// Flag that indicates if the overall convex hull needs to be recomputed (i.e. mass trace convex hulls were modified)
    mutable bool convex_hulls_modified_{};

    /// Overall convex hull of the feature
    mutable ConvexHull2D convex_hull_;

    /// subordinate features (e.g. features that represent alternative explanations, usually with lower quality)
    std::vector<Feature> subordinates_;

  };

} // namespace OpenMS
