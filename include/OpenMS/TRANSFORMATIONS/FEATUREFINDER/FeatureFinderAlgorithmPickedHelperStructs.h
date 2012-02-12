// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

#include <vector>

namespace OpenMS
{

  /**
   * @brief Wrapper struct for all the classes needed by the FeatureFinderAlgorithmPicked and the associated classes
   *
   * @see FeatureFinderAlgorithmPicked
   * @see TraceFitter
   */
  struct OPENMS_DLLAPI FeatureFinderAlgorithmPickedHelperStructs
  {

    /**
     * @brief Helper structure for seeds used in FeatureFinderAlgorithmPicked
     */
    struct OPENMS_DLLAPI Seed
    {
      ///Spectrum index
      Size spectrum;
      ///Peak index
      Size peak;
      ///Intensity
      Real intensity;

      /// Comparison operator
      bool operator<(const Seed& rhs) const
      {
        return intensity<rhs.intensity;
      }
    };

    /**
     * @brief Helper struct for mass traces used in FeatureFinderAlgorithmPicked
     */
    template<class PeakType>
    struct MassTrace
    {
      ///Maximum peak pointer
      const PeakType* max_peak;
      ///RT of maximum peak
      DoubleReal max_rt;

      ///Theoretical intensity value (scaled to [0,1])
      DoubleReal theoretical_int;

      ///Contained peaks (pair of RT and pointer to peak)
      std::vector<std::pair<DoubleReal, const PeakType*> > peaks;

      ///determines the convex hull of the trace
      ConvexHull2D getConvexhull() const
      {
        ConvexHull2D::PointArrayType hull_points(peaks.size());
        for (Size i=0; i<peaks.size(); ++i)
        {
          hull_points[i][0] = peaks[i].first;
          hull_points[i][1] = peaks[i].second->getMZ();
        }
				ConvexHull2D hull;
				hull.addPoints(hull_points);
        return hull;
      }

      ///Sets the maximum to the highest contained peak of the trace
      void updateMaximum()
      {
        if (peaks.empty()) return;

        max_rt = peaks.begin()->first;
        max_peak = peaks.begin()->second;

        for (Size i=1; i<peaks.size(); ++i)
        {
          if (peaks[i].second->getIntensity()>max_peak->getIntensity())
          {
            max_rt = peaks[i].first;
            max_peak = peaks[i].second;
          }
        }
      }

      ///Returns the average m/z of all peaks in this trace (weighted by intensity)
      DoubleReal getAvgMZ() const
      {
        DoubleReal sum = 0.0;
        DoubleReal intensities = 0.0;
        for (Size i=0; i<peaks.size(); ++i)
        {
          sum += peaks[i].second->getMZ()*peaks[i].second->getIntensity();
          intensities += peaks[i].second->getIntensity();
        }
        return sum / intensities;
      }

      ///Checks if this Trace is valid (has more than 2 points)
      bool isValid() const
      {
        return (peaks.size()>=3);
      }

    };

    /**
     * @brief Helper struct for a collection of mass traces used in FeatureFinderAlgorithmPicked
     */
    template<class PeakType>
    struct MassTraces
      : public std::vector< MassTrace<PeakType> >
    {
      /// Constructor
      MassTraces()
        : max_trace(0)
      {
      }

      /// Returns the peak count of all traces
      Size getPeakCount() const
      {
        Size sum = 0;
        for (Size i=0; i<this->size(); ++i)
        {
          sum += this->at(i).peaks.size();
        }
        return sum;
      }

      ///Checks if still valid (seed still contained and enough traces)
      bool isValid(DoubleReal seed_mz, DoubleReal trace_tolerance)
      {
        //Abort if too few traces were found
        if (this->size()<2) return false;

        //Abort if the seed was removed
        for (Size j=0; j<this->size(); ++j)
        {
          if (std::fabs(seed_mz-this->at(j).getAvgMZ())<=trace_tolerance)
          {
            return true;
          }
        }
        return false;
      }

      /**
        @brief Returns the theoretical maximum trace index

        @exception Exception::Precondition is thrown if there are not mass traces (not only in debug mode)
      */
      Size getTheoreticalmaxPosition() const
      {
        if (!this->size())
        {
          throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There must be at least one trace to determine the theoretical maximum trace!");
        }

        Size max=0;
        DoubleReal max_int=this->at(0).theoretical_int;
        for (Size i=1; i<this->size(); ++i)
        {
          if (this->at(i).theoretical_int>max_int)
          {
            max_int = this->at(i).theoretical_int;
            max = i;
          }
        }
        return max;
      }

      ///Sets the baseline to the lowest contained peak of the trace
      void updateBaseline()
      {
        if (this->size() == 0)
        {
          baseline = 0.0;
          return;
        }
        bool first = true;
        for (Size i=0 ; i < this->size() ; ++i)
        {
          for (Size j=0 ; j < this->at(i).peaks.size() ; ++j)
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

      /**
        @brief Returns the RT boundaries of the mass traces

        @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
      */
      std::pair<DoubleReal,DoubleReal> getRTBounds() const
      {
        if (!this->size())
        {
          throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There must be at least one trace to determine the RT boundaries!");
        }

        DoubleReal min = std::numeric_limits<DoubleReal>::max();
        DoubleReal max = -std::numeric_limits<DoubleReal>::max();
        //Abort if the seed was removed
        for (Size i=0; i<this->size(); ++i)
        {
          for (Size j=0; j<this->at(i).peaks.size(); ++j)
          {
            DoubleReal rt = this->at(i).peaks[j].first;
            if (rt>max) max = rt;
            if (rt<min) min = rt;
          }
        }
        return std::make_pair(min,max);
      }

      /// Maximum intensity trace
      Size max_trace;
      /// Estimated baseline in the region of the feature (used for the fit)
      DoubleReal baseline;
    };

    /**
     * @brief Helper structure for a theoretical isotope pattern used in FeatureFinderAlgorithmPicked
     */
    struct OPENMS_DLLAPI TheoreticalIsotopePattern
    {
      ///Vector of intensity contributions
      std::vector<DoubleReal> intensity;
      ///Number of optional peaks at the beginning of the pattern
      Size optional_begin;
      ///Number of optional peaks at the end of the pattern
      Size optional_end;
      ///The maximum intensity contribution before scaling the pattern to 1
      DoubleReal max;
      ///The number of isotopes trimmed on the left side. This is needed to reconstruct the monoisotopic peak.
      Size trimmed_left;
      /// Returns the size
      Size size() const
      {
        return intensity.size();
      }
    };

    /**
     * @brief Helper structure for a found isotope pattern used in FeatureFinderAlgorithmPicked
     */
    struct OPENMS_DLLAPI IsotopePattern
    {
      ///Peak index (-1 if peak was not found, -2 if it was removed to improve the isotope fit)
      std::vector<SignedSize> peak;
      ///Spectrum index (undefined if peak index is -1 or -2)
      std::vector<Size> spectrum;
      ///Peak intensity (0 if peak index is -1 or -2)
      std::vector<DoubleReal> intensity;
      ///m/z score of peak (0 if peak index is -1 or -2)
      std::vector<DoubleReal> mz_score;
      ///Theoretical m/z value of the isotope peak
      std::vector<DoubleReal> theoretical_mz;
      ///Theoretical isotope pattern
      TheoreticalIsotopePattern theoretical_pattern;

      /// Constructor that resizes the internal vectors
      IsotopePattern(Size size)
        : peak(size,-1),
          spectrum(size),
          intensity(size),
          mz_score(size),
          theoretical_mz(size)
      {
      }
    };

  };
}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H
