// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  FeatureFileOptions::FeatureFileOptions()
    : loadConvexhull_(true),
      loadSubordinates_(true),
      metadata_only_(false),
      has_rt_range_(false),
      has_mz_range_(false),
      has_intensity_range_(false)
  {
  }
  
  FeatureFileOptions::~FeatureFileOptions()
  {
  }
  
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

  void FeatureFileOptions::setRTRange(const DRange<1>& range)
  {
    rt_range_ = range;
    has_rt_range_ = true;
  }
  
  bool FeatureFileOptions::hasRTRange() const
  {
    return has_rt_range_;
  }
  
  const DRange<1>& FeatureFileOptions::getRTRange() const
  {
    return rt_range_;
  }
  
  void FeatureFileOptions::setMZRange(const DRange<1>& range)
  {
    mz_range_ = range;
    has_mz_range_ = true;
  }

  bool FeatureFileOptions::hasMZRange() const
  {
    return has_mz_range_;
  }
  
  const DRange<1>& FeatureFileOptions::getMZRange() const
  {
    return mz_range_;
  }
  
  void FeatureFileOptions::setIntensityRange(const DRange<1>& range)
  {
    intensity_range_ = range;
    has_intensity_range_ = true;
  }
  
  bool FeatureFileOptions::hasIntensityRange() const
  {
    return has_intensity_range_;
  }
  
  const DRange<1>& FeatureFileOptions::getIntensityRange() const
  {
    return intensity_range_;
  }
  
} // namespace OpenMS
