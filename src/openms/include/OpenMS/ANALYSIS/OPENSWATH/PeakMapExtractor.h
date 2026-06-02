// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <limits>
#include <vector>

namespace OpenMS
{
  /**
    @brief Extract raw mz/RT/IM peak clouds for targeted OpenSWATH coordinates.

    The extractor follows the same targeting contract as
    @ref ChromatogramExtractorAlgorithm: it consumes sorted extraction
    coordinates and applies RT, m/z, and optionally ion mobility windows.
    Unlike chromatogram extraction, it retains every matching raw point instead
    of integrating a scalar intensity per spectrum.

    The output is one extracted peak map per coordinate. Each peak map stores
    parallel arrays for m/z, RT, ion mobility, and intensity.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI PeakMapExtractor :
    public ProgressLogger
  {
  public:
    using ExtractionCoordinates = ChromatogramExtractorAlgorithm::ExtractionCoordinates;

    /**
      @brief One extracted targeted peak map.
    */
    struct OPENMS_DLLAPI ExtractedPeakMap
    {
      String native_id;
      double target_mz{0.0};
      double target_rt{std::numeric_limits<double>::quiet_NaN()};
      double target_ion_mobility{-1.0};
      double rt_start{0.0};
      double rt_end{-1.0};

      std::vector<double> mz;
      std::vector<double> rt;
      std::vector<double> ion_mobility;
      std::vector<double> intensity;
    };

    /**
      @brief Extract raw peak maps at the coordinates defined by @p extraction_coordinates.

      @param[in] input Input spectral map
      @param[out] output Output peak maps; one row per extraction coordinate
      @param[in] extraction_coordinates Target coordinates to extract
      @param[in] mz_extraction_window Full m/z extraction window in Th or ppm
      @param[in] ppm Whether @p mz_extraction_window is specified in ppm
      @param[in] im_extraction_window Full ion mobility extraction window. A
        positive value activates ion mobility filtering around the coordinate's
        target ion mobility. A negative value extracts the full ion mobility
        range while still recording the raw ion mobility values.
      @param[in] filter Extraction profile in m/z space. Only @c "tophat" is
        supported because the output preserves raw points instead of applying
        weighted integration.
    */
    void extractPeakMaps(const OpenSwath::SpectrumAccessPtr& input,
                         std::vector<ExtractedPeakMap>& output,
                         const std::vector<ExtractionCoordinates>& extraction_coordinates,
                         double mz_extraction_window,
                         bool ppm,
                         double im_extraction_window,
                         const String& filter = String("tophat"));

  private:
    int getFilterNr_(const String& filter) const;
  };
} // namespace OpenMS
