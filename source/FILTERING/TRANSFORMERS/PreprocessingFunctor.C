// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
