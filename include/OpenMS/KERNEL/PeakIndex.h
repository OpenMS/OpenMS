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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_PEAKINDEX_H
#define OPENMS_KERNEL_PEAKINDEX_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <limits>

namespace OpenMS
{
  /**
    @brief Index of a peak or feature

    This struct can be used to store both peak or feature indices.
  */
  struct PeakIndex
  {
    /// Default constructor. Creates an invalid peak reference
    inline PeakIndex() :
      peak((std::numeric_limits<Size>::max)()),
      spectrum((std::numeric_limits<Size>::max)())
    {}

    /// Constructor that sets the peak index (for feaure maps)
    explicit inline PeakIndex(Size peak) :
      peak(peak),
      spectrum((std::numeric_limits<Size>::max)())
    {}

    /// Constructor that sets the peak and spectrum index (for peak maps)
    inline PeakIndex(Size spectrum, Size peak) :
      peak(peak),
      spectrum(spectrum)
    {}

    /// returns if the current peak ref is valid
    inline bool isValid() const
    {
      return peak != (std::numeric_limits<Size>::max)();
    }

    /// Invalidates the current index
    inline void clear()
    {
      peak = (std::numeric_limits<Size>::max)();
      spectrum = (std::numeric_limits<Size>::max)();
    }

    /**
      @brief Access to the feature (or consensus feature) corresponding to this index

      This method is intended for arrays of features e.g. FeatureMap

      The main advantage of using this method instead accessing the data directly is that range
      check performed in debug mode.

      @exception Exception::Precondition is thrown if this index is invalid for the @p map (only in
      debug mode)
    */
    template <typename FeatureMapType>
    const typename FeatureMapType::value_type & getFeature(const FeatureMapType & map) const
    {
      OPENMS_PRECONDITION(peak < map.size(), "Feature index exceeds map size");
      return map[peak];
    }

    /**
      @brief Access to a peak corresponding to this index.

      This method is intended for arrays of DSpectra e.g. MSExperiment

      The main advantage of using this method instead accessing the data directly is that range
      check performed in debug mode.

      @exception Exception::Precondition is thrown if this index is invalid for the @p map (only in
      debug mode)
    */
    template <typename PeakMapType>
    const typename PeakMapType::PeakType & getPeak(const PeakMapType & map) const
    {
      OPENMS_PRECONDITION(spectrum < map.size(), "Spectrum index exceeds map size");
      OPENMS_PRECONDITION(peak < map[spectrum].size(), "Peak index exceeds spectrum size");
      return map[spectrum][peak];
    }

    /**
      @brief Access to a spectrum corresponding to this index

      This method is intended for arrays of DSpectra e.g. MSExperiment.

      The main advantage of using this method instead accessing the data directly is that range
      check performed in debug mode.

      @exception Exception::Precondition is thrown if this index is invalid for the @p map (only in
      debug mode)
    */
    template <typename PeakMapType>
    const typename PeakMapType::SpectrumType & getSpectrum(const PeakMapType & map) const
    {
      OPENMS_PRECONDITION(spectrum < map.size(), "Spectrum index exceeds map size");
      return map[spectrum];
    }

    /// Equality operator
    inline bool operator==(const PeakIndex & rhs) const
    {
      return peak == rhs.peak && spectrum == rhs.spectrum;
    }

    /// Inequality operator
    inline bool operator!=(const PeakIndex & rhs) const
    {
      return peak != rhs.peak || spectrum != rhs.spectrum;
    }

    /// Peak or feature index
    Size peak;
    /// Spectrum index
    Size spectrum;
  };

} // namespace OpenMS

#endif // OPENMS_KERNEL_PEAKINDEX_H
