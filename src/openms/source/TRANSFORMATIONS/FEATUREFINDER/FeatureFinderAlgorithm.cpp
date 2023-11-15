// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMRM.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  void FeatureFinderAlgorithm::registerChildren()
  {
    Factory<FeatureFinderAlgorithm>::registerProduct
    (
      FeatureFinderAlgorithmPicked::getProductName(),
      &FeatureFinderAlgorithmPicked::create
    );
    Factory<FeatureFinderAlgorithm>::registerProduct
    (
      FeatureFinderAlgorithmMRM::getProductName(),
      &FeatureFinderAlgorithmMRM::create
    );

  }

}
