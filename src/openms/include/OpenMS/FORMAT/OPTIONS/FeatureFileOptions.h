// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Options for loading files containing features.

    A passive value object: each setter stores its value, each getter
    returns the stored value. The defaults loaded by the constructor are
    @c loadConvexHull == @c true, @c loadSubordinates == @c true,
    @c metadataOnly == @c false, @c sizeOnly == @c false, and no
    RT / m/z / intensity range set (i.e. the @c has*Range queries return
    @c false).

    Reader implementations interpret these options. The @c has*Range
    flags indicate whether a range was ever supplied --
    @c getRTRange / @c getMZRange / @c getIntensityRange only return a
    meaningful range when the corresponding @c has*Range is @c true.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureFileOptions
  {
public:
    /// Default constructor; see class brief for the default option values.
    FeatureFileOptions();
    /// Destructor.
    ~FeatureFileOptions();

    ///@name Convex-hull option
    ///@{
    /// Set whether to load each feature's convex hull (default @c true).
    void setLoadConvexHull(bool convex);
    /// Whether to load each feature's convex hull.
    bool getLoadConvexHull() const;
    ///@}

    ///@name Subordinate option
    ///@{
    /// Set whether to load a feature's subordinates (default @c true).
    void setLoadSubordinates(bool sub);
    /// Whether to load a feature's subordinates.
    bool getLoadSubordinates() const;
    ///@}

    ///@name Metadata-only option
    ///@{
    /// Set whether to skip the per-feature payload and load only metadata (default @c false).
    void setMetadataOnly(bool only);
    /// Whether only metadata is loaded.
    bool getMetadataOnly() const;
    ///@}

    ///@name Size-only option
    ///@{
    /// Set whether to read just the feature count and skip the actual features (default @c false).
    void setSizeOnly(bool only);
    /// Whether to read just the feature count.
    bool getSizeOnly() const;
    ///@}

    ///@name RT-range option
    ///@{
    /// Restrict loaded features to the given RT range. Also marks the range as set (@ref hasRTRange returns @c true afterwards).
    void setRTRange(const DRange<1> & range);
    /// @c true if @ref setRTRange has been called.
    bool hasRTRange() const;
    /// The configured RT range; only meaningful when @ref hasRTRange is @c true.
    const DRange<1> & getRTRange() const;
    ///@}

    ///@name m/z range option
    ///@{
    /// Restrict loaded features to the given m/z range. Also marks the range as set (@ref hasMZRange returns @c true afterwards).
    void setMZRange(const DRange<1> & range);
    /// @c true if @ref setMZRange has been called.
    bool hasMZRange() const;
    /// The configured m/z range; only meaningful when @ref hasMZRange is @c true.
    const DRange<1> & getMZRange() const;
    ///@}

    ///@name Intensity-range option
    ///@{
    /// Restrict loaded features to the given intensity range. Also marks the range as set (@ref hasIntensityRange returns @c true afterwards).
    void setIntensityRange(const DRange<1> & range);
    /// @c true if @ref setIntensityRange has been called.
    bool hasIntensityRange() const;
    /// The configured intensity range; only meaningful when @ref hasIntensityRange is @c true.
    const DRange<1> & getIntensityRange() const;
    ///@}

private:
    bool loadConvexhull_;
    bool loadSubordinates_;
    bool metadata_only_;
    bool has_rt_range_;
    bool has_mz_range_;
    bool has_intensity_range_;
    bool size_only_;
    DRange<1> rt_range_;
    DRange<1> mz_range_;
    DRange<1> intensity_range_;

  };

} // namespace OpenMS

