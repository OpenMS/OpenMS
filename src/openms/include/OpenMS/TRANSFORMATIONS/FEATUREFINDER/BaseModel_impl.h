// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <OpenMS/SIMULATION/EGHModel.h>

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
    Factory<BaseModel<1> >::registerProduct(EGHModel::getProductName(), &EGHModel::create);

    return;
  }

} // namespace OpenMS

