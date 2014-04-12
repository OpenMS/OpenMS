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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/MATH/MISC/Spline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	SplinePackage::SplinePackage(std::vector<double> mz, std::vector<double> intensity, double scaling) : spline_(3, mz, intensity)
	{				
        if (!(mz.size() == intensity.size() && mz.size() > 2))
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"m/z and intensity vectors either not of the same size or too short.");
        }
        	
		mzMin_ = *min_element(mz.begin(), mz.end());
		mzMax_ = *max_element(mz.begin(), mz.end());
		mzStepWidth_ = scaling*(mzMax_ - mzMin_)/(mz.size() - 1);    // step width somewhat smaller than the average raw data spacing	
	}
	
	SplinePackage::~SplinePackage() 
	{
	}
	
	double SplinePackage::getMzMin()
	{
		return mzMin_;
	}

	double SplinePackage::getMzMax()
	{
		return mzMax_;
	}
	
	double SplinePackage::getMzStepWidth()
	{
		return mzStepWidth_;
	}
	
	bool SplinePackage::isInPackage(double mz)
	{
		return (mz >= mzMin_ && mz <= mzMax_);
	}
	
	Spline2d<double> SplinePackage::getSpline()
	{
		return spline_;
	}
	
	double SplinePackage::eval(double mz)
	{
		if (this->isInPackage(mz))
		{
			double intensity = spline_.eval(mz);
			if (intensity < 0)
			{
				return 0;
			}
			else
			{
				return intensity;
			}
		}
		else
		{
			return 0;
		}
	}

}
