// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
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
  PreprocessingFunctor::PreprocessingFunctor() :
    DefaultParamHandler("PreprocessingFunctor")
  {
  }

  PreprocessingFunctor::PreprocessingFunctor(const PreprocessingFunctor & source) :
    DefaultParamHandler(source)
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

  PreprocessingFunctor & PreprocessingFunctor::operator=(const PreprocessingFunctor & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
