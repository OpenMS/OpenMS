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
  (for peptides). The spatial extent of the Feature in the LC-MS data set is
  represented by a bounding box (RT and m/z range). Optionally, a list of
  convex hulls (one for each isotopic peak) provides per-trace detail. The
  model description can store the parameters of a two-dimensional theoretical
  model of the underlying signal in LC-MS. Currently, non-peptide compounds
  are also represented as features.

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

    ///@name Bounding box
    //@{
    /// Returns the bounding box of the feature in RT (dim 0) and m/z (dim 1)
    const DBoundingBox<2>& getBoundingBox() const;
    /// Sets the bounding box of the feature
    void setBoundingBox(const DBoundingBox<2>& box);
    /// Sets the bounding box from RT and m/z ranges
    void setBoundingBox(double rt_min, double mz_min, double rt_max, double mz_max);
    /// Returns true if a bounding box has been set
    bool hasBoundingBox() const;
    /// Compute and set the bounding box from the convex hulls; returns true if updated
    bool updateBoundingBoxFromConvexHulls();
    //@}

    ///@name Mass traces (convex hulls)
    //@{
    /// Non-mutable access to the convex hulls
    const std::vector<ConvexHull2D>& getConvexHulls() const;
    /// Set the convex hulls of single mass traces
    void setConvexHulls(const std::vector<ConvexHull2D>& hulls);
    /// Add a single convex hull
    void addConvexHull(const ConvexHull2D& hull);
    /// Remove all convex hulls
    void clearConvexHulls();

    /// Returns if the mass trace convex hulls (or bounding box) enclose the position specified by @p rt and @p mz
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

    /// Bounding box in RT (dim 0) and m/z (dim 1)
    DBoundingBox<2> bounding_box_;

    /// Array of convex hulls (one for each mass trace, optional)
    std::vector<ConvexHull2D> convex_hulls_;

    /// subordinate features (e.g. features that represent alternative explanations, usually with lower quality)
    std::vector<Feature> subordinates_;

  };

} // namespace OpenMS
