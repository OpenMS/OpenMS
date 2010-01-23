// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_IMPL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_IMPL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>

// include derived classes here
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaIsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/SIMULATION/EGHModel.h>

namespace OpenMS
{
  
    template<>
    void BaseModel<2>::registerChildren()
    {
        Factory< BaseModel<2> >::registerProduct(ProductModel<2>::getProductName(), &ProductModel<2>::create);
    }
  
    template<>
    void BaseModel<1>::registerChildren()
    {

        Factory< BaseModel<1> >::registerProduct(GaussModel::getProductName(), &GaussModel::create);
        Factory< BaseModel<1> >::registerProduct(BiGaussModel::getProductName(), &BiGaussModel::create);
        Factory< BaseModel<1> >::registerProduct(IsotopeModel::getProductName(), &IsotopeModel::create);
        Factory< BaseModel<1> >::registerProduct(ExtendedIsotopeModel::getProductName(), &ExtendedIsotopeModel::create);
        Factory< BaseModel<1> >::registerProduct(LmaIsotopeModel::getProductName(), &LmaIsotopeModel::create);
        Factory< BaseModel<1> >::registerProduct(EmgModel::getProductName(), &EmgModel::create);
        Factory< BaseModel<1> >::registerProduct(LmaGaussModel::getProductName(), &LmaGaussModel::create);
        Factory< BaseModel<1> >::registerProduct(EGHModel::getProductName(), &EGHModel::create);

        return;
    }

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_IMPL_H

