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
  */
  class OPENMS_DLLAPI FeatureFileOptions
  {
public:
    ///Default constructor
    FeatureFileOptions();
    ///Destructor
    ~FeatureFileOptions();

    ///@name convex hull option
    ///sets whether or not to load convex hull
    void setLoadConvexHull(bool convex);
    ///returns whether or not to load convex hull
    bool getLoadConvexHull() const;

    ///@name subordinate option
    ///sets whether or not load subordinates
    void setLoadSubordinates(bool sub);
    ///returns whether or not to load subordinates
    bool getLoadSubordinates() const;

    ///@name metadata option
    ///sets whether or not to load only meta data
    void setMetadataOnly(bool only);
    ///returns whether or not to load only meta data
    bool getMetadataOnly() const;

    ///@name lazyload option
    ///sets whether or not to load only feature count
    void setSizeOnly(bool only);
    ///returns whether or not to load only meta data
    bool getSizeOnly() const;

    ///@name RT range option
    ///restricts the range of RT values for peaks to load
    void setRTRange(const DRange<1> & range);
    ///returns @c true if an RT range has been set
    bool hasRTRange() const;
    ///returns the RT range
    const DRange<1> & getRTRange() const;

    ///@name m/z range option
    ///restricts the range of MZ values for peaks to load
    void setMZRange(const DRange<1> & range);
    ///returns @c true if an MZ range has been set
    bool hasMZRange() const;
    ///returns the MZ range
    const DRange<1> & getMZRange() const;

    ///@name Intensity range option
    ///restricts the range of intensity values for peaks to load
    void setIntensityRange(const DRange<1> & range);
    ///returns @c true if an intensity range has been set
    bool hasIntensityRange() const;
    ///returns the intensity range
    const DRange<1> & getIntensityRange() const;

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

