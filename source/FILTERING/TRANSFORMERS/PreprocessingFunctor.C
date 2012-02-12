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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>
#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  PreprocessingFunctor::PreprocessingFunctor()
    : DefaultParamHandler("PreprocessingFunctor")
  {
  }

  PreprocessingFunctor::PreprocessingFunctor(const PreprocessingFunctor& source)
    : DefaultParamHandler(source)
  {
		
  }

	PreprocessingFunctor::~PreprocessingFunctor()
	{
	}
	
	void PreprocessingFunctor::registerChildren()
	{
    Factory<PreprocessingFunctor>::registerProduct(ThresholdMower::getProductName(), &ThresholdMower::create);
    Factory<PreprocessingFunctor>::registerProduct(WindowMower::getProductName(), &WindowMower::create);
    Factory<PreprocessingFunctor>::registerProduct(Scaler::getProductName(), &Scaler::create);
    Factory<PreprocessingFunctor>::registerProduct(NLargest::getProductName(), &NLargest::create);
    Factory<PreprocessingFunctor>::registerProduct(MarkerMower::getProductName(), &MarkerMower::create);
    Factory<PreprocessingFunctor>::registerProduct(SqrtMower::getProductName(), &SqrtMower::create);
    Factory<PreprocessingFunctor>::registerProduct(Normalizer::getProductName(), &Normalizer::create);
    Factory<PreprocessingFunctor>::registerProduct(ParentPeakMower::getProductName(), &ParentPeakMower::create);
		Factory<PreprocessingFunctor>::registerProduct(BernNorm::getProductName(), &BernNorm::create);
	}
	
  PreprocessingFunctor& PreprocessingFunctor::operator = (const PreprocessingFunctor& source)
  {
		if (this != &source)
		{
   		DefaultParamHandler::operator=(source);
		}
    return *this;
  }
}
