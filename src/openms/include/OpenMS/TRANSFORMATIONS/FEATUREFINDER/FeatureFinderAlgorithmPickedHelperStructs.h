// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Stephan Aiche $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>

#include <vector>
#include <list>
#include <cmath>

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
      float intensity;

      /// Comparison operator
      bool operator<(const Seed& rhs) const;

    };

    /**
     * @brief Helper struct for mass traces used in FeatureFinderAlgorithmPicked
     */
    template <class PeakType>
    struct MassTrace
    {
      ///Maximum peak pointer
      const PeakType* max_peak;
      ///RT of maximum peak
      double max_rt;

      ///Theoretical intensity value (scaled to [0,1])
      double theoretical_int;

      ///Contained peaks (pair of RT and pointer to peak)
      std::vector<std::pair<double, const PeakType*> > peaks;

      ///determines the convex hull of the trace
      ConvexHull2D getConvexhull() const
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

      ///Sets the maximum to the highest contained peak of the trace
      void updateMaximum()
      {
        if (peaks.empty()) return;

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

      ///Returns the average m/z of all peaks in this trace (weighted by intensity)
      double getAvgMZ() const
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

      ///Checks if this Trace is valid (has more than 2 points)
      bool isValid() const
      {
        return peaks.size() >= 3;
      }

    };

    /**
     * @brief Helper struct for a collection of mass traces used in FeatureFinderAlgorithmPicked
     */
    template <class PeakType>
    struct MassTraces :
      private std::vector<MassTrace<PeakType> >
    {
      typedef std::vector<MassTrace<PeakType> > privvec;

      // public exports of used methods
      using privvec::size;
      using privvec::at;
      using privvec::reserve;
      using privvec::push_back;
      using privvec::operator[];
      using privvec::back;
      using privvec::clear;
      using privvec::begin;
      using privvec::end;
      typedef typename privvec::iterator iterator;
      typedef typename privvec::const_iterator const_iterator;

      /// Constructor
      MassTraces() :
        max_trace(0)
      {
      }

      /// Returns the peak count of all traces
      Size getPeakCount() const
      {
        Size sum = 0;
        for (Size i = 0; i < this->size(); ++i)
        {
          sum += this->at(i).peaks.size();
        }
        return sum;
      }

      ///Checks if still valid (seed still contained and enough traces)
      bool isValid(double seed_mz, double trace_tolerance)
      {
        //Abort if too few traces were found
        if (this->size() < 2) return false;

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

      /**
        @brief Returns the theoretical maximum trace index

        @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
      */
      Size getTheoreticalmaxPosition() const
      {
        if (!this->size())
        {
          throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There must be at least one trace to determine the theoretical maximum trace!");
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

      ///Sets the baseline to the lowest contained peak of the trace
      void updateBaseline()
      {
        if (this->size() == 0)
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

      /**
        @brief Returns the RT boundaries of the mass traces

        @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
      */
      std::pair<double, double> getRTBounds() const
      {
        if (!this->size())
        {
          throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There must be at least one trace to determine the RT boundaries!");
        }

        double min = std::numeric_limits<double>::max();
        double max = -std::numeric_limits<double>::max();
        //Abort if the seed was removed
        for (Size i = 0; i < this->size(); ++i)
        {
          for (Size j = 0; j < this->at(i).peaks.size(); ++j)
          {
            double rt = this->at(i).peaks[j].first;
            if (rt > max) max = rt;
            if (rt < min) min = rt;
          }
        }
        return std::make_pair(min, max);
      }

      /**
        @brief Computes a flat representation of MassTraces, i.e., a single
               intensity value for each point in RT. The flattened representation
               is comparable to the TIC of the MassTraces.

        @param intensity_profile An empty std::list of pair<double, double> that will be filled.
                The first element of the pair holds the RT value, the second value the sum of intensities
                of all peaks in the different mass traces with this specific RT.
      */
      void computeIntensityProfile(std::list<std::pair<double, double> >& intensity_profile) const
      {
        // typedefs for better readability
        typedef typename MassTraces<PeakType>::const_iterator TTraceIterator;
        typedef std::list<std::pair<double, double> >::iterator TProfileIterator;
        typedef typename std::vector<std::pair<double, const PeakType*> > TMassTracePeakList;
        typedef typename TMassTracePeakList::const_iterator TTracePeakIterator;

        TTraceIterator trace_it = this->begin();
        // we add the first trace without check, as the profile is currently empty
        for (TTracePeakIterator trace_peak_it = trace_it->peaks.begin(); trace_peak_it != trace_it->peaks.end(); ++trace_peak_it)
        {
          intensity_profile.push_back(std::make_pair(trace_peak_it->first, trace_peak_it->second->getIntensity()));
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
              intensity_profile.push_back(std::make_pair(trace_peak_it->first, trace_peak_it->second->getIntensity()));
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

      /// Maximum intensity trace
      Size max_trace;
      /// Estimated baseline in the region of the feature (used for the fit)
      double baseline;
    };

    /**
     * @brief Helper structure for a theoretical isotope pattern used in FeatureFinderAlgorithmPicked
     */
    struct OPENMS_DLLAPI TheoreticalIsotopePattern
    {
      ///Vector of intensity contributions
      std::vector<double> intensity;
      ///Number of optional peaks at the beginning of the pattern
      Size optional_begin;
      ///Number of optional peaks at the end of the pattern
      Size optional_end;
      ///The maximum intensity contribution before scaling the pattern to 1
      double max;
      ///The number of isotopes trimmed on the left side. This is needed to reconstruct the monoisotopic peak.
      Size trimmed_left;
      /// Returns the size
      Size size() const;

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
      std::vector<double> intensity;
      ///m/z score of peak (0 if peak index is -1 or -2)
      std::vector<double> mz_score;
      ///Theoretical m/z value of the isotope peak
      std::vector<double> theoretical_mz;
      ///Theoretical isotope pattern
      TheoreticalIsotopePattern theoretical_pattern;

      /// Constructor that resizes the internal vectors
      explicit IsotopePattern(Size size);

    };

  };
}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKEDHELPERSTRUCTS_H
