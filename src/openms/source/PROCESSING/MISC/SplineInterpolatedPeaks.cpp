// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>
#include <OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  SplineInterpolatedPeaks::SplineInterpolatedPeaks(const std::vector<double>& pos, const std::vector<double>& intensity)
  {
    SplineInterpolatedPeaks::init_(pos, intensity);
  }

  SplineInterpolatedPeaks::SplineInterpolatedPeaks(const MSSpectrum& raw_spectrum)
  {
    std::vector<double> pos;
    std::vector<double> intensity;
    for (const auto &it : raw_spectrum)
    {
      pos.push_back(it.getMZ());
      intensity.push_back(it.getIntensity());
    }
    SplineInterpolatedPeaks::init_(pos, intensity);
  }

  SplineInterpolatedPeaks::SplineInterpolatedPeaks(const MSChromatogram& raw_chromatogram)
  {
    std::vector<double> rt;
    std::vector<double> intensity;
    for (const auto &it : raw_chromatogram)
    {
      rt.push_back(it.getRT());
      intensity.push_back(it.getIntensity());
    }
    SplineInterpolatedPeaks::init_(rt, intensity);
  }

  SplineInterpolatedPeaks::~SplineInterpolatedPeaks() = default;

  void SplineInterpolatedPeaks::init_(const std::vector<double>& pos, const std::vector<double>& intensity)
  {

    if (!(pos.size() == intensity.size() && pos.size() > 2))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z and intensity vectors either not of the same size or too short.");
    }

    const double new_package = 2; // start a new package if delta m/z is greater than new_package times previous one

    pos_min_ = pos.front();
    pos_max_ = pos.back();

    // remove unnecessary zeros, i.e. zero intensity data points with zeros to the left and right
    std::vector<double> pos_slim1; // slimmer vector after removal of zero-intensity datapoints from pos
    std::vector<double> intensity_slim1; // slimmer vector after removal of zero-intensity datapoints from intensity
    pos_slim1.reserve(pos.size());
    intensity_slim1.reserve(intensity.size());
    if (intensity[0] != 0 || intensity[1] != 0)
    {
      pos_slim1.push_back(pos[0]);
      intensity_slim1.push_back(intensity[0]);
    }
    bool last_intensity_zero = (intensity[0] == 0);
    bool current_intensity_zero = (intensity[0] == 0);
    bool next_intensity_zero = (intensity[1] == 0);
    for (size_t i = 1; i < pos.size() - 1; ++i)
    {
      last_intensity_zero = current_intensity_zero;
      current_intensity_zero = next_intensity_zero;
      next_intensity_zero = (intensity[i + 1] == 0);
      if (!last_intensity_zero || !current_intensity_zero || !next_intensity_zero)
      {
        pos_slim1.push_back(pos[i]);
        intensity_slim1.push_back(intensity[i]);
      }
    }
    if (intensity[pos.size() - 1] != 0 || intensity[pos.size() - 2] != 0)
    {
      pos_slim1.push_back(pos[pos.size() - 1]);
      intensity_slim1.push_back(intensity[pos.size() - 1]);
    }

    // remove Thermo bug zeros
    // (In some Thermo data appear odd zero intensity data points. Normal data points are sometimes quickly followed by a zero.
    // These zeros are clearly not part of the profile, but bugs. The following code snippet removes them. A datapoint is
    // "quickly followed" by a second one, if the m/z step is shorter than scaling_Thermo_bug times the previous m/z step.)
    std::vector<double> pos_slim2; // slimmer vector after removal of Thermo bugs from pos_slim1
    std::vector<double> intensity_slim2; // slimmer vector after removal of Thermo bugs from intensity_slim1
    const double scaling_Thermo_bug = 1.0 / 50.0; // scaling factor for Thermo bug
    pos_slim2.reserve(pos_slim1.size());
    intensity_slim2.reserve(intensity_slim1.size());
    pos_slim2.push_back(pos_slim1[0]);
    pos_slim2.push_back(pos_slim1[1]);
    intensity_slim2.push_back(intensity_slim1[0]);
    intensity_slim2.push_back(intensity_slim1[1]);
    for (size_t i = 2; i < pos_slim1.size(); ++i)
    {
      if (intensity_slim1[i] == 0)
      {
        if ((pos_slim1[i] - pos_slim1[i - 1]) < (pos_slim1[i - 1] - pos_slim1[i - 2]) * scaling_Thermo_bug)
        {
          continue;
        }
      }
      pos_slim2.push_back(pos_slim1[i]);
      intensity_slim2.push_back(intensity_slim1[i]);
    }

    // subdivide spectrum into packages
    std::vector<bool> start_package;
    start_package.push_back(true);
    start_package.push_back(false);
    for (size_t i = 2; i < pos_slim2.size(); ++i)
    {
      start_package.push_back((pos_slim2[i] - pos_slim2[i - 1]) / (pos_slim2[i - 1] - pos_slim2[i - 2]) > new_package);
    }

    // fill the packages
    std::vector<double> pos_package;
    std::vector<double> intensity_package;
    for (size_t i = 0; i < pos_slim2.size(); ++i)
    {
      if (start_package[i] && i > 0)
      {
        if (intensity_package.size() > 1)
        {
          // Two or more data points in package. At least one of them will be non-zero since unnecessary zeros removed above.
          packages_.emplace_back(pos_package, intensity_package);
        }
        pos_package.clear();
        intensity_package.clear();
      }
      pos_package.push_back(pos_slim2[i]);
      intensity_package.push_back(intensity_slim2[i]);
    }
    // add the last package
    if (intensity_package.size() > 1)
    {
      packages_.emplace_back(pos_package, intensity_package);
    }

  }

  double SplineInterpolatedPeaks::getPosMin() const
  {
    return pos_min_;
  }

  double SplineInterpolatedPeaks::getPosMax() const
  {
    return pos_max_;
  }

  size_t SplineInterpolatedPeaks::size() const
  {
    return packages_.size();
  }

  SplineInterpolatedPeaks::Navigator SplineInterpolatedPeaks::getNavigator(double scaling)
  {
    if (packages_.empty())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0, "packages_ must not be empty");
    }
    return Navigator(&packages_, pos_max_, scaling);
  }

  SplineInterpolatedPeaks::Navigator::Navigator(const std::vector<SplinePackage>* packages, double pos_max, double scaling) :
    packages_(packages),
    last_package_(0),
    pos_max_(pos_max),
    pos_step_width_scaling_(scaling)
  {
  }

  SplineInterpolatedPeaks::Navigator::Navigator() = default;

  SplineInterpolatedPeaks::Navigator::~Navigator() = default;

  size_t SplineInterpolatedPeaks::Navigator::findPackage_(double pos) const
  {
    // The packages are ordered and separated by gaps (init_() splits sorted positions into them), so the package
    // containing pos is the last one that starts at or before pos.
    const std::vector<SplinePackage>& packages = *packages_;
    const SplinePackage& last = packages[last_package_];
    if (pos >= last.getPosMin())
    {
      // the last package, or one of the next two (positions are usually requested in increasing order)
      if (pos <= last.getPosMax() || last_package_ + 1 == packages.size() || pos < packages[last_package_ + 1].getPosMin())
      {
        return last_package_;
      }
      if (last_package_ + 2 == packages.size() || pos < packages[last_package_ + 2].getPosMin())
      {
        return last_package_ + 1;
      }
    }
    const auto right = std::upper_bound(packages.begin(), packages.end(), pos,
                                        [](double p, const SplinePackage& package) { return p < package.getPosMin(); });
    return right == packages.begin() ? packages.size() : static_cast<size_t>(right - packages.begin()) - 1;
  }

  double SplineInterpolatedPeaks::Navigator::eval(double pos)
  {
    const size_t i = findPackage_(pos);
    if (i == packages_->size())
    {
      // left of the first package
      last_package_ = 0;
      return 0.0;
    }
    last_package_ = i;
    // zero in the gap after package i (eval() checks the range)
    return (*packages_)[i].eval(pos);
  }

  double SplineInterpolatedPeaks::Navigator::getNextPos(double pos)
  {
    const std::vector<SplinePackage>& packages = *packages_;
    const size_t max_index = packages.size() - 1;
    const size_t i = findPackage_(pos);

    if (i == packages.size())
    {
      // left of the first package: jump to its start
      last_package_ = 0;
      return packages[0].getPosMin();
    }
    const SplinePackage& package = packages[i];
    if (pos > package.getPosMax())
    {
      // in the gap after package i: jump to the start of the next package (or to the end)
      if (i == max_index)
      {
        last_package_ = max_index;
        return pos_max_;
      }
      last_package_ = i + 1;
      return packages[i + 1].getPosMin();
    }

    // find m/z in the package
    if (pos + pos_step_width_scaling_ * package.getPosStepWidth() > package.getPosMax())
    {
      // The next step gets us outside the current package.
      // Let's move to the package to the right.
      // check index limit
      if (i == max_index)
      {
        last_package_ = max_index;
        return pos_max_;
      }
      // jump to min m/z of next package
      last_package_ = i + 1;
      return packages[i + 1].getPosMin();
    }
    else
    {
      // make a small step within the package
      last_package_ = i;
      return pos + pos_step_width_scaling_ * package.getPosStepWidth();
    }
  }

}
