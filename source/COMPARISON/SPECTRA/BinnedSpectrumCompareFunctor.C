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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/CONCEPT/Factory.h>



using namespace std;

namespace OpenMS
{
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor()
		: DefaultParamHandler(BinnedSpectrumCompareFunctor::getProductName())
  {
  }
  
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source)
		:	DefaultParamHandler(source)
  {
  }
 
	BinnedSpectrumCompareFunctor::~BinnedSpectrumCompareFunctor()
	{
	}
 
  BinnedSpectrumCompareFunctor& BinnedSpectrumCompareFunctor::operator=(const BinnedSpectrumCompareFunctor& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }
 
  void BinnedSpectrumCompareFunctor::registerChildren()
  {
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSharedPeakCount::getProductName(), &BinnedSharedPeakCount::create);
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSpectralContrastAngle::getProductName(), &BinnedSpectralContrastAngle::create);
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSumAgreeingIntensities::getProductName(), &BinnedSumAgreeingIntensities::create);
  }
  
  BinnedSpectrumCompareFunctor::IncompatibleBinning::IncompatibleBinning(const char* file, int line, const char* function, const char* message) throw()
          : BaseException(file, line, function, "BinnedSpectrumCompareFunctor::IncompatibleBinning", message)
  {
  } 
  BinnedSpectrumCompareFunctor::IncompatibleBinning::~IncompatibleBinning() throw()
  {	
  }

}
