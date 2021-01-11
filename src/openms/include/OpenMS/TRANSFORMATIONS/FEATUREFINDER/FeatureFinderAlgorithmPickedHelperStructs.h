// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/KERNEL/Peak1D.h>

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
    struct OPENMS_DLLAPI MassTrace
    {
      ///Maximum peak pointer
      const Peak1D* max_peak;
      ///RT of maximum peak
      double max_rt;

      ///Theoretical intensity value (scaled to [0,1])
      double theoretical_int;

      ///Contained peaks (pair of RT and pointer to peak)
      std::vector<std::pair<double, const Peak1D*> > peaks;

      ///determines the convex hull of the trace
      ConvexHull2D getConvexhull() const;

      ///Sets the maximum to the highest contained peak of the trace
      void updateMaximum();

      ///Returns the average m/z of all peaks in this trace (weighted by intensity)
      double getAvgMZ() const;

      ///Checks if this Trace is valid (has more than 2 points)
      bool isValid() const;

    };

    /**
     * @brief Helper struct for a collection of mass traces used in FeatureFinderAlgorithmPicked
     */
    struct OPENMS_DLLAPI MassTraces :
      private std::vector<MassTrace>
    {
      typedef std::vector<MassTrace> privvec;

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
      typedef privvec::iterator iterator;
      typedef privvec::const_iterator const_iterator;

      /// Constructor
      MassTraces();

      /// Returns the peak count of all traces
      Size getPeakCount() const;

      ///Checks if still valid (seed still contained and enough traces)
      bool isValid(double seed_mz, double trace_tolerance);

      /**
        @brief Returns the theoretical maximum trace index

        @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
      */
      Size getTheoreticalmaxPosition() const;

      ///Sets the baseline to the lowest contained peak of the trace
      void updateBaseline();

      /**
        @brief Returns the RT boundaries of the mass traces

        @exception Exception::Precondition is thrown if there are no mass traces (not only in debug mode)
      */
      std::pair<double, double> getRTBounds() const;

      /**
        @brief Computes a flat representation of MassTraces, i.e., a single
               intensity value for each point in RT. The flattened representation
               is comparable to the TIC of the MassTraces.

        @param intensity_profile An empty std::list of pair<double, double> that will be filled.
                The first element of the pair holds the RT value, the second value the sum of intensities
                of all peaks in the different mass traces with this specific RT.
      */
      void computeIntensityProfile(std::list<std::pair<double, double> >& intensity_profile) const;

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

