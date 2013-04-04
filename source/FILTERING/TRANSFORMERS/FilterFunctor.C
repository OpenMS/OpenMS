// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  FilterFunctor::FilterFunctor() :
    DefaultParamHandler("FilterFunctor")
  {
  }

  FilterFunctor::FilterFunctor(const FilterFunctor & source) :
    DefaultParamHandler(source)
  {
  }

  FilterFunctor & FilterFunctor::operator=(const FilterFunctor & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  FilterFunctor::~FilterFunctor()
  {
  }

  void FilterFunctor::registerChildren()
  {
    Factory<FilterFunctor>::registerProduct(ComplementFilter::getProductName(), &ComplementFilter::create);
    Factory<FilterFunctor>::registerProduct(GoodDiffFilter::getProductName(), &GoodDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(IntensityBalanceFilter::getProductName(), &IntensityBalanceFilter::create);
    Factory<FilterFunctor>::registerProduct(NeutralLossDiffFilter::getProductName(), &NeutralLossDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(IsotopeDiffFilter::getProductName(), &IsotopeDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(TICFilter::getProductName(), &TICFilter::create);
  }

}
