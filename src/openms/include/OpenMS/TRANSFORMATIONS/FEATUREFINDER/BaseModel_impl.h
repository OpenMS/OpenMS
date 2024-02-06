// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>

// include derived classes here
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>

namespace OpenMS
{
  
  OPENMS_DLLAPI void BaseModel::registerChildren()
  {
    Factory<BaseModel>::registerProduct(GaussModel::getProductName(), &GaussModel::create);
    Factory<BaseModel>::registerProduct(BiGaussModel::getProductName(), &BiGaussModel::create);
    Factory<BaseModel>::registerProduct(IsotopeModel::getProductName(), &IsotopeModel::create);
    Factory<BaseModel>::registerProduct(ExtendedIsotopeModel::getProductName(), &ExtendedIsotopeModel::create);
    Factory<BaseModel>::registerProduct(EmgModel::getProductName(), &EmgModel::create);

    return;
  }

} // namespace OpenMS

