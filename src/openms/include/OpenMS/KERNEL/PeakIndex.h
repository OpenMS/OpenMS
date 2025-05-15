// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once


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

    /// Constructor that sets the peak index (for feature maps)
    explicit inline PeakIndex(Size lpeak) :
      peak(lpeak),
      spectrum((std::numeric_limits<Size>::max)())
    {}

    /// Constructor that sets the peak and spectrum index (for peak maps)
    inline PeakIndex(Size lspectrum, Size lpeak) :
      peak(lpeak),
      spectrum(lspectrum)
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

