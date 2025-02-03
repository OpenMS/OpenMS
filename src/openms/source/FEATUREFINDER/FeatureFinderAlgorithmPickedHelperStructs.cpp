// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

namespace OpenMS
{
  FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern::IsotopePattern(Size size) :
    peak(size, -1),
    spectrum(size),
    intensity(size),
    mz_score(size),
    theoretical_mz(size)
  {
  }

  Size FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern::size() const
  {
    return intensity.size();
  }

  bool FeatureFinderAlgorithmPickedHelperStructs::Seed::operator<(const Seed& rhs) const
  {
    return intensity < rhs.intensity;
  }

  ConvexHull2D FeatureFinderAlgorithmPickedHelperStructs::MassTrace::getConvexhull() const
  {
    ConvexHull2D::PointArrayType hull_points(peaks.size());
    for (Size i = 0; i < peaks.size(); ++i)
    {
      hull_points[i][0] = peaks[i].first;
      hull_points[i][1] = peaks[i].second->getMZ();
    }
    ConvexHull2D hull;
    hull.addPoints(hull_points);
    return hull;
  }

  void FeatureFinderAlgorithmPickedHelperStructs::MassTrace::updateMaximum()
  {
    if (peaks.empty())
    {
      return;
    }
    max_rt = peaks.begin()->first;
    max_peak = peaks.begin()->second;

    for (Size i = 1; i < peaks.size(); ++i)
    {
      if (peaks[i].second->getIntensity() > max_peak->getIntensity())
      {
        max_rt = peaks[i].first;
        max_peak = peaks[i].second;
      }
    }
  }

  double FeatureFinderAlgorithmPickedHelperStructs::MassTrace::getAvgMZ() const
  {
    double sum = 0.0;
    double intensities = 0.0;
    for (Size i = 0; i < peaks.size(); ++i)
    {
      sum += peaks[i].second->getMZ() * peaks[i].second->getIntensity();
      intensities += peaks[i].second->getIntensity();
    }
    return sum / intensities;
  }

  bool FeatureFinderAlgorithmPickedHelperStructs::MassTrace::isValid() const
  {
    return peaks.size() >= 3;
  }

  FeatureFinderAlgorithmPickedHelperStructs::MassTraces::MassTraces() :
    max_trace(0)
  {
  }

  Size FeatureFinderAlgorithmPickedHelperStructs::MassTraces::getPeakCount() const
  {
    Size sum = 0;
    for (Size i = 0; i < this->size(); ++i)
    {
      sum += this->at(i).peaks.size();
    }
    return sum;
  }

  bool FeatureFinderAlgorithmPickedHelperStructs::MassTraces::isValid(double seed_mz, double trace_tolerance)
  {
    //Abort if too few traces were found
    if (this->size() < 2)
    {
      return false;
    }
    //Abort if the seed was removed
    for (Size j = 0; j < this->size(); ++j)
    {
      if (std::fabs(seed_mz - this->at(j).getAvgMZ()) <= trace_tolerance)
      {
        return true;
      }
    }
    return false;
  }

  Size FeatureFinderAlgorithmPickedHelperStructs::MassTraces::getTheoreticalmaxPosition() const
  {
    if (this->empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one trace to determine the theoretical maximum trace!");
    }

    Size max = 0;
    double max_int = this->at(0).theoretical_int;
    for (Size i = 1; i < this->size(); ++i)
    {
      if (this->at(i).theoretical_int > max_int)
      {
        max_int = this->at(i).theoretical_int;
        max = i;
      }
    }
    return max;
  }

  void FeatureFinderAlgorithmPickedHelperStructs::MassTraces::updateBaseline()
  {
    if (this->empty())
    {
      baseline = 0.0;
      return;
    }
    bool first = true;
    for (Size i = 0; i < this->size(); ++i)
    {
      for (Size j = 0; j < this->at(i).peaks.size(); ++j)
      {
        if (first)
        {
          baseline = this->at(i).peaks[j].second->getIntensity();
          first = false;
        }
        if (this->at(i).peaks[j].second->getIntensity() < baseline)
        {
          baseline = this->at(i).peaks[j].second->getIntensity();
        }
      }
    }
  }

  std::pair<double, double> FeatureFinderAlgorithmPickedHelperStructs::MassTraces::getRTBounds() const
  {
    if (this->empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one trace to determine the RT boundaries!");
    }

    double min = std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::max();
    //Abort if the seed was removed
    for (Size i = 0; i < this->size(); ++i)
    {
      for (Size j = 0; j < this->at(i).peaks.size(); ++j)
      {
        double rt = this->at(i).peaks[j].first;
        if (rt > max)
        {
          max = rt;
        }
        if (rt < min)
        {
          min = rt;
        }
      }
    }
    return std::make_pair(min, max);
  }

  void FeatureFinderAlgorithmPickedHelperStructs::MassTraces::computeIntensityProfile(std::list<std::pair<double, double> >& intensity_profile) const
  {
    // typedefs for better readability
    typedef MassTraces::const_iterator TTraceIterator;
    typedef std::list<std::pair<double, double> >::iterator TProfileIterator;
    typedef std::vector<std::pair<double, const Peak1D*> > TMassTracePeakList;
    typedef TMassTracePeakList::const_iterator TTracePeakIterator;

    TTraceIterator trace_it = this->begin();
    // we add the first trace without check, as the profile is currently empty
    for (TTracePeakIterator trace_peak_it = trace_it->peaks.begin(); trace_peak_it != trace_it->peaks.end(); ++trace_peak_it)
    {
      intensity_profile.emplace_back(trace_peak_it->first, trace_peak_it->second->getIntensity());
    }
    ++trace_it;

    // accumulate intensities over all the remaining mass traces
    for (; trace_it != this->end(); ++trace_it)
    {
      TProfileIterator profile_it = intensity_profile.begin();
      TTracePeakIterator trace_peak_it = trace_it->peaks.begin();

      while (trace_peak_it != trace_it->peaks.end())
      {
        // append .. if profile has already ended
        if (profile_it == intensity_profile.end())
        {
          intensity_profile.emplace_back(trace_peak_it->first, trace_peak_it->second->getIntensity());
          ++trace_peak_it;
        }
        // prepend
        else if (profile_it->first > trace_peak_it->first)
        {
          intensity_profile.insert(profile_it, std::make_pair(trace_peak_it->first, trace_peak_it->second->getIntensity()));
          ++trace_peak_it;
        }
        // proceed
        else if (profile_it->first < trace_peak_it->first)
        {
          ++profile_it;
        }
        // merge
        else if (profile_it->first == trace_peak_it->first)
        {
          profile_it->second += trace_peak_it->second->getIntensity();
          ++trace_peak_it;
          ++profile_it;
        }
      }
    }
  }

}
