// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>

// include derived classes here
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>

namespace OpenMS
{

  template <>
  OPENMS_DLLAPI void BaseModel<2>::registerChildren()
  {
    Factory<BaseModel<2> >::registerProduct(ProductModel<2>::getProductName(), &ProductModel<2>::create);
  }

  template <>
  OPENMS_DLLAPI void BaseModel<1>::registerChildren()
  {

    Factory<BaseModel<1> >::registerProduct(GaussModel::getProductName(), &GaussModel::create);
    Factory<BaseModel<1> >::registerProduct(BiGaussModel::getProductName(), &BiGaussModel::create);
    Factory<BaseModel<1> >::registerProduct(IsotopeModel::getProductName(), &IsotopeModel::create);
    Factory<BaseModel<1> >::registerProduct(ExtendedIsotopeModel::getProductName(), &ExtendedIsotopeModel::create);
    Factory<BaseModel<1> >::registerProduct(EmgModel::getProductName(), &EmgModel::create);

    return;
  }

} // namespace OpenMS

