// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>

using namespace std;

namespace OpenMS
{
  FeatureFileOptions::FeatureFileOptions() :
    loadConvexhull_(true),
    loadSubordinates_(true),
    metadata_only_(false),
    has_rt_range_(false),
    has_mz_range_(false),
    has_intensity_range_(false),
    size_only_(false)
  {
  }

  FeatureFileOptions::~FeatureFileOptions() = default;

  void FeatureFileOptions::setLoadConvexHull(bool convex)
  {
    loadConvexhull_ = convex;
  }

  bool FeatureFileOptions::getLoadConvexHull() const
  {
    return loadConvexhull_;
  }

  void FeatureFileOptions::setLoadSubordinates(bool sub)
  {
    loadSubordinates_ = sub;
  }

  bool FeatureFileOptions::getLoadSubordinates() const
  {
    return loadSubordinates_;
  }

  void FeatureFileOptions::setMetadataOnly(bool only)
  {
    metadata_only_ = only;
  }

  bool FeatureFileOptions::getMetadataOnly() const
  {
    return metadata_only_;
  }

  void FeatureFileOptions::setRTRange(const DRange<1> & range)
  {
    rt_range_ = range;
    has_rt_range_ = true;
  }

  bool FeatureFileOptions::hasRTRange() const
  {
    return has_rt_range_;
  }

  const DRange<1> & FeatureFileOptions::getRTRange() const
  {
    return rt_range_;
  }

  void FeatureFileOptions::setMZRange(const DRange<1> & range)
  {
    mz_range_ = range;
    has_mz_range_ = true;
  }

  bool FeatureFileOptions::hasMZRange() const
  {
    return has_mz_range_;
  }
  
  bool FeatureFileOptions::getSizeOnly() const
  {
    return size_only_;
  }

  void FeatureFileOptions::setSizeOnly(bool size_only)
  {
    size_only_ = size_only;
  }
  
  const DRange<1> & FeatureFileOptions::getMZRange() const
  {
    return mz_range_;
  }

  void FeatureFileOptions::setIntensityRange(const DRange<1> & range)
  {
    intensity_range_ = range;
    has_intensity_range_ = true;
  }

  bool FeatureFileOptions::hasIntensityRange() const
  {
    return has_intensity_range_;
  }

  const DRange<1> & FeatureFileOptions::getIntensityRange() const
  {
    return intensity_range_;
  }

} // namespace OpenMS
