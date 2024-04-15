// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFinder.h>

#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  FeatureFinder::FeatureFinder() :
    flags_()
  {
  }

  FeatureFinder::~FeatureFinder() = default;

  Param FeatureFinder::getParameters(const String& /*algorithm_name*/) const
  {/* TODO: remove?
    Param tmp;
    if (algorithm_name != "none")
    {
      FeatureFinderAlgorithm* a = Factory<FeatureFinderAlgorithm>::create(algorithm_name);
      tmp.insert("", a->getDefaultParameters());
      delete(a);
    }
    return tmp;
    */
   std::cerr << "TODO: remove" << std::endl;
   return Param();
  }

}
