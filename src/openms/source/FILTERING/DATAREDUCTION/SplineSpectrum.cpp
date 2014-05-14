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
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/MATH/MISC/Spline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

SplineSpectrum::SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity)
{
	SplineSpectrum::init_(mz, intensity, 0.7);
}

SplineSpectrum::SplineSpectrum(const std::vector<double>& mz, const std::vector<double>& intensity, double scaling)
{
	SplineSpectrum::init_(mz, intensity, scaling);
}

SplineSpectrum::SplineSpectrum(MSSpectrum<Peak1D> rawSpectrum)
{
	std::vector<double> mz;
	std::vector<double> intensity;
	for (MSSpectrum<Peak1D>::Iterator it = rawSpectrum.begin(); it != rawSpectrum.end(); ++it)
	{
		mz.push_back(it->getMZ());
		intensity.push_back(it->getIntensity());
	}
	SplineSpectrum::init_(mz, intensity, 0.7);
}

SplineSpectrum::SplineSpectrum(MSSpectrum<Peak1D> rawSpectrum, double scaling)
{
	std::vector<double> mz;
	std::vector<double> intensity;
	for (MSSpectrum<Peak1D>::Iterator it = rawSpectrum.begin(); it != rawSpectrum.end(); ++it)
	{
		mz.push_back(it->getMZ());
		intensity.push_back(it->getIntensity());
	}
	SplineSpectrum::init_(mz, intensity, scaling);
}

SplineSpectrum::~SplineSpectrum()
{
}

void SplineSpectrum::init_(const std::vector<double>& mz, const std::vector<double>& intensity, double scaling)
{

	if (!(mz.size() == intensity.size() && mz.size() > 2))
	{
		throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"m/z and intensity vectors either not of the same size or too short.");
	}

	const double newPackage = 2;    // start a new package if delta m/z is greater than newPackage times previous one

	mzMin_ = mz[0];
	mzMax_ = mz[mz.size() -1];

	// remove unnecessary zeros, i.e. zero intensity data points with zeros to the left and right
	std::vector<double> mzSlim;
	std::vector<double> intensitySlim;
	if (intensity[0]!=0 || intensity[1]!=0)
	{
		mzSlim.push_back(mz[0]);
		intensitySlim.push_back(intensity[0]);
	}
	bool lastIntensityZero = (intensity[0] == 0);
	bool currentIntensityZero = (intensity[0] == 0);
	bool nextIntensityZero = (intensity[1] == 0);
	for (unsigned i=1; i<mz.size()-1; ++i)
	{
		lastIntensityZero = currentIntensityZero;
		currentIntensityZero = nextIntensityZero;
		nextIntensityZero = (intensity[i+1] == 0);
		if (!lastIntensityZero || !currentIntensityZero || !nextIntensityZero)
		{
			mzSlim.push_back(mz[i]);
			intensitySlim.push_back(intensity[i]);
		}
	}
	if (intensity[mz.size()-1]!=0 || intensity[mz.size()-2]!=0)
	{
		mzSlim.push_back(mz[mz.size()-1]);
		intensitySlim.push_back(intensity[mz.size()-1]);
	}

	// subdivide spectrum into packages
	std::vector<bool> startPackage;
	startPackage.push_back(true);
	startPackage.push_back(false);
	for (unsigned i=2; i<mzSlim.size(); ++i)
	{
		startPackage.push_back((mzSlim[i] - mzSlim[i-1])/(mzSlim[i-1] - mzSlim[i-2]) > newPackage);
	}

	// fill the packages
	std::vector<double> mzPackage;
	std::vector<double> intensityPackage;
	for (unsigned i=0; i<mzSlim.size(); ++i)
	{
		if (startPackage[i] && i > 0)
		{
			if (intensityPackage.size() > 2)
			{
				// Three or more data points in package. At least one of them will be non-zero since unnecessary zeros removed above.
				packages_.push_back(SplinePackage(mzPackage, intensityPackage, scaling));
			}
			mzPackage.clear();
			intensityPackage.clear();
		}
		mzPackage.push_back(mzSlim[i]);
		intensityPackage.push_back(intensitySlim[i]);
	}
	// add the last package
	if (intensityPackage.size() > 2)
	{
		packages_.push_back(SplinePackage(mzPackage, intensityPackage, scaling));
	}
}

double SplineSpectrum::getMzMin() const
{
	return mzMin_;
}

double SplineSpectrum::getMzMax() const
{
	return mzMax_;
}

SplinePackage SplineSpectrum::getPackage(int i) const
{
	return packages_[i];
}

SplineSpectrum::Navigator::Navigator(const std::vector<SplinePackage> * packages, double mzMin, double mzMax) : packages_(packages), lastPackage_(0), mzMin_(mzMin), mzMax_(mzMax)
{
}

SplineSpectrum::Navigator::~Navigator()
{
}

double SplineSpectrum::Navigator::eval(double mz)
{
	if (mz < (*packages_)[lastPackage_].getMzMin())
	{
		for (int i = lastPackage_; i >= 0; --i)
		{
			if (mz > (*packages_)[i].getMzMax())
			{
				lastPackage_ = i;
				return 0.0;
			}
			if (mz >= (*packages_)[i].getMzMin())
			{
				lastPackage_ = i;
				return (*packages_)[i].eval(mz);
			}
		}
	}
	else
	{
		for (int i = lastPackage_; i < (int)(*packages_).size(); ++i)
		{
			if (mz < (*packages_)[i].getMzMin())
			{
				lastPackage_ = i;
				return 0.0;
			}
			if (mz <= (*packages_)[i].getMzMax())
			{
				lastPackage_ = i;
				return (*packages_)[i].eval(mz);
			}
		}
	}
	return 0.0;
}

double SplineSpectrum::Navigator::getNextMz(double mz)
{

	int minIndex = 0;
	int maxIndex = (*packages_).size() - 1;
	int i = lastPackage_;
	SplinePackage package = (*packages_)[i];

	// find correct package
	while (!(package.isInPackage(mz)))
	{
		if (mz < package.getMzMin())
		{
			--i;
			// check index limit
			if (i < minIndex)
			{
				lastPackage_ = minIndex;
				return (*packages_)[minIndex].getMzMin();
			}
			// m/z in the gap?
			package = (*packages_)[i];
			if (mz > package.getMzMax())
			{
				lastPackage_ = i + 1;
				return (*packages_)[i+1].getMzMin();
			}
		}
		else if (mz > package.getMzMax())
		{

			++i;
			// check index limit
			if (i > maxIndex)
			{
				lastPackage_ = maxIndex;
				return mzMax_;
			}
			// m/z in the gap?
			package = (*packages_)[i];
			if (mz < package.getMzMin())
			{
				lastPackage_ = i;
				return package.getMzMin();
			}
		}
	}

	// find m/z in the package
	if (mz + package.getMzStepWidth() > package.getMzMax())
	{
		// The next step gets us outside the current package.
		// Let's move to the package to the right.
		++i;
		// check index limit
		if (i > maxIndex)
		{
			lastPackage_ = maxIndex;
			return mzMax_;
		}
		// jump to min m/z of next package
		lastPackage_ = i;
		return (*packages_)[i].getMzMin();
	}
	else
	{
		// make a small step within the package
		lastPackage_ = i;
		return mz + package.getMzStepWidth();
	}
}

SplineSpectrum::Navigator SplineSpectrum::getNavigator()
{
	return Navigator(&packages_, mzMin_, mzMax_);
}

}
