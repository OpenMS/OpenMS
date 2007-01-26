// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <iostream>

// all from BaseModel derived classes
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>




namespace OpenMS
{

	template<>
	void BaseModel<2>::registerChildren(){
		Factory< BaseModel<2> >::registerProduct(ProductModel<2>::getProductName(), &ProductModel<2>::create);
	}

	template<>
	void BaseModel<1>::registerChildren(){

		Factory< BaseModel<1> >::registerProduct(GaussModel::getProductName(), &GaussModel::create);
		Factory< BaseModel<1> >::registerProduct(BiGaussModel::getProductName(), &BiGaussModel::create);
		Factory< BaseModel<1> >::registerProduct(IsotopeModel::getProductName(), &IsotopeModel::create);
		Factory< BaseModel<1> >::registerProduct(EmgModel::getProductName(), &EmgModel::create);
		Factory< BaseModel<1> >::registerProduct(LmaGaussModel::getProductName(), &LmaGaussModel::create);
		Factory< BaseModel<1> >::registerProduct(LogNormalModel::getProductName(), &LogNormalModel::create);

		return;
	}

} // namespace OpenMS
