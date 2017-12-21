// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

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

  SplineSpectrum::SplineSpectrum(MSSpectrum& raw_spectrum)
  {
    std::vector<double> mz;
    std::vector<double> intensity;
    for (MSSpectrum::Iterator it = raw_spectrum.begin(); it != raw_spectrum.end(); ++it)
    {
      mz.push_back(it->getMZ());
      intensity.push_back(it->getIntensity());
    }
    SplineSpectrum::init_(mz, intensity, 0.7);
  }

  SplineSpectrum::SplineSpectrum(MSSpectrum& raw_spectrum, double scaling)
  {
    std::vector<double> mz;
    std::vector<double> intensity;
    for (MSSpectrum::Iterator it = raw_spectrum.begin(); it != raw_spectrum.end(); ++it)
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z and intensity vectors either not of the same size or too short.");
    }

    const double new_package = 2; // start a new package if delta m/z is greater than new_package times previous one

    mz_min_ = mz.front();
    mz_max_ = mz.back();

    // remove unnecessary zeros, i.e. zero intensity data points with zeros to the left and right
    std::vector<double> mz_slim1; // slimmer vector after removal of zero-intensity datapoints from mz
    std::vector<double> intensity_slim1; // slimmer vector after removal of zero-intensity datapoints from intensity
    mz_slim1.reserve(mz.size());
    intensity_slim1.reserve(intensity.size());
    if (intensity[0] != 0 || intensity[1] != 0)
    {
      mz_slim1.push_back(mz[0]);
      intensity_slim1.push_back(intensity[0]);
    }
    bool last_intensity_zero = (intensity[0] == 0);
    bool current_intensity_zero = (intensity[0] == 0);
    bool next_intensity_zero = (intensity[1] == 0);
    for (size_t i = 1; i < mz.size() - 1; ++i)
    {
      last_intensity_zero = current_intensity_zero;
      current_intensity_zero = next_intensity_zero;
      next_intensity_zero = (intensity[i + 1] == 0);
      if (!last_intensity_zero || !current_intensity_zero || !next_intensity_zero)
      {
        mz_slim1.push_back(mz[i]);
        intensity_slim1.push_back(intensity[i]);
      }
    }
    if (intensity[mz.size() - 1] != 0 || intensity[mz.size() - 2] != 0)
    {
      mz_slim1.push_back(mz[mz.size() - 1]);
      intensity_slim1.push_back(intensity[mz.size() - 1]);
    }

    // remove Thermo bug zeros
    // (In some Thermo data appear odd zero intensity data points. Normal data points are sometimes quickly followed by a zero.
    // These zeros are clearly not part of the profile, but bugs. The following code snippet removes them. A datapoint is
    // "quickly followed" by a second one, if the m/z step is shorter than scaling_Thermo_bug times the previous m/z step.)
    std::vector<double> mz_slim2; // slimmer vector after removal of Thermo bugs from mz_slim1
    std::vector<double> intensity_slim2; // slimmer vector after removal of Thermo bugs from intensity_slim1
    const double scaling_Thermo_bug = 1.0 / 50.0; // scaling factor for Thermo bug
    mz_slim2.reserve(mz_slim1.size());
    intensity_slim2.reserve(intensity_slim1.size());
    mz_slim2.push_back(mz_slim1[0]);
    mz_slim2.push_back(mz_slim1[1]);
    intensity_slim2.push_back(intensity_slim1[0]);
    intensity_slim2.push_back(intensity_slim1[1]);
    for (size_t i = 2; i < mz_slim1.size(); ++i)
    {
      if (intensity_slim1[i] == 0)
      {
        if ((mz_slim1[i] - mz_slim1[i - 1]) < (mz_slim1[i - 1] - mz_slim1[i - 2]) * scaling_Thermo_bug)
        {
          continue;
        }
      }
      mz_slim2.push_back(mz_slim1[i]);
      intensity_slim2.push_back(intensity_slim1[i]);
    }

    // subdivide spectrum into packages
    std::vector<bool> start_package;
    start_package.push_back(true);
    start_package.push_back(false);
    for (size_t i = 2; i < mz_slim2.size(); ++i)
    {
      start_package.push_back((mz_slim2[i] - mz_slim2[i - 1]) / (mz_slim2[i - 1] - mz_slim2[i - 2]) > new_package);
    }

    // fill the packages
    std::vector<double> mz_package;
    std::vector<double> intensity_package;
    for (size_t i = 0; i < mz_slim2.size(); ++i)
    {
      if (start_package[i] && i > 0)
      {
        if (intensity_package.size() > 1)
        {
          // Two or more data points in package. At least one of them will be non-zero since unnecessary zeros removed above.
          packages_.push_back(SplinePackage(mz_package, intensity_package, scaling));
        }
        mz_package.clear();
        intensity_package.clear();
      }
      mz_package.push_back(mz_slim2[i]);
      intensity_package.push_back(intensity_slim2[i]);
    }
    // add the last package
    if (intensity_package.size() > 1)
    {
      packages_.push_back(SplinePackage(mz_package, intensity_package, scaling));
    }

  }

  double SplineSpectrum::getMzMin() const
  {
    return mz_min_;
  }

  double SplineSpectrum::getMzMax() const
  {
    return mz_max_;
  }

  size_t SplineSpectrum::getSplineCount() const
  {
    return packages_.size();
  }

  SplineSpectrum::Navigator SplineSpectrum::getNavigator()
  {
    if (packages_.empty())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
    }
    return Navigator(&packages_, mz_min_, mz_max_);
  }

  SplineSpectrum::Navigator::Navigator(const std::vector<SplinePackage>* packages, double mz_min, double mz_max) :
    packages_(packages),
    last_package_(0),
    mz_min_(mz_min),
    mz_max_(mz_max)
  {
  }

  SplineSpectrum::Navigator::Navigator()
  {
  }

  SplineSpectrum::Navigator::~Navigator()
  {
  }

  double SplineSpectrum::Navigator::eval(double mz)
  {
    if (mz < (*packages_)[last_package_].getMzMin())
    { // look left
      for (int i = (int) last_package_; i >= 0; --i)
      {
        if (mz > (*packages_)[i].getMzMax())
        {
          last_package_ = i;
          return 0.0;
        }
        if (mz >= (*packages_)[i].getMzMin())
        {
          last_package_ = i;
          return (*packages_)[i].eval(mz);
        }
      }
    }
    else
    { // look right
      for (size_t i = last_package_; i < (size_t)(*packages_).size(); ++i)
      {
        if (mz < (*packages_)[i].getMzMin())
        {
          last_package_ = i;
          return 0.0;
        }
        if (mz <= (*packages_)[i].getMzMax())
        {
          last_package_ = i;
          return (*packages_)[i].eval(mz);
        }
      }
    }
    return 0.0;
  }

  double SplineSpectrum::Navigator::getNextMz(double mz)
  {

    int min_index = 0;
    int max_index = static_cast<Int>((*packages_).size()) - 1;
    int i = static_cast<Int>(last_package_);
    SplinePackage package = (*packages_)[i];

    // find correct package
    while (!(package.isInPackage(mz)))
    {
      if (mz < package.getMzMin())
      {
        --i;
        // check index limit
        if (i < min_index)
        {
          last_package_ = min_index;
          return (*packages_)[min_index].getMzMin();
        }
        // m/z in the gap?
        package = (*packages_)[i];
        if (mz > package.getMzMax())
        {
          last_package_ = i + 1;
          return (*packages_)[i + 1].getMzMin();
        }
      }
      else if (mz > package.getMzMax())
      {

        ++i;
        // check index limit
        if (i > max_index)
        {
          last_package_ = max_index;
          return mz_max_;
        }
        // m/z in the gap?
        package = (*packages_)[i];
        if (mz < package.getMzMin())
        {
          last_package_ = i;
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
      if (i > max_index)
      {
        last_package_ = max_index;
        return mz_max_;
      }
      // jump to min m/z of next package
      last_package_ = i;
      return (*packages_)[i].getMzMin();
    }
    else
    {
      // make a small step within the package
      last_package_ = i;
      return mz + package.getMzStepWidth();
    }
  }

}
