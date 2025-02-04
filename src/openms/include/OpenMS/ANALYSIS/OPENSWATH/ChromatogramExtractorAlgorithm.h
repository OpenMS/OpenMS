// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

namespace OpenMS
{

  /**
   * @brief The ChromatogramExtractorAlgorithm extracts chromatograms from a MS data.
   *
   * It will take as input a set of transitions coordinates and will extract
   * the signal of the provided map at the product ion m/z and retention time
   * (rt) values specified by the extraction coordinates. This interface only
   * expects a set of coordinates which are up to the user to fill but a
   * convenient prepare_coordinates function is provided (in the
   * ChromatogramExtractor class) to create the coordinates for the most common
   * case of an MS2 and MS1 extraction.
   *
   * In the case of MS2 extraction, the map is assumed to originate from a SWATH
   * (data-independent acquisition or DIA) experiment.
   *
  */
  class OPENMS_DLLAPI ChromatogramExtractorAlgorithm :
    public ProgressLogger
  {

public:

    struct ExtractionCoordinates
    {
      double mz = 0.0; ///< m/z value around which should be extracted
      double ion_mobility = 0.0; ///< ion mobility value around which should be extracted
      double mz_precursor = 0.0; ///< precursor m/z value (is currently ignored by the algorithm)
      double rt_start = 0.0; ///< rt start of extraction (in seconds)
      double rt_end = 0.0; ///< rt end of extraction (in seconds)
      std::string id; ///< identifier

      static bool SortExtractionCoordinatesByMZ(
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& left,
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& right)
      {
        return left.mz < right.mz;
      }
      static bool SortExtractionCoordinatesReverseByMZ(
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& left,
          const ChromatogramExtractorAlgorithm::ExtractionCoordinates& right)
      {
        return left.mz > right.mz;
      }
    };

    /**
     * @brief Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates.
     *
     * @param input Input spectral map
     * @param output Output chromatograms (XICs)
     * @param extraction_coordinates Extracts around these coordinates (from
     *   rt_start to rt_end in seconds - extracts the whole chromatogram if
     *   rt_end - rt_start < 0).
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension in Th or ppm (e.g. a window of 50 ppm means an extraction of
     * 25 ppm on either side)
     * @param ppm Whether mz_extraction_window is in ppm or in Th
     * @param im_extraction_window Full window width (i.e. twice the tolerance) for IM extraction. Must be positive.
     * @param filter Which function to apply in m/z space (currently "tophat" only)
     *
    */
    void extractChromatograms(const OpenSwath::SpectrumAccessPtr& input,
        std::vector< OpenSwath::ChromatogramPtr >& output,
        const std::vector<ExtractionCoordinates>& extraction_coordinates,
        double mz_extraction_window,
        bool ppm,
        double im_extraction_window,
        const String& filter);

    /**
     * @brief Extract the next mz value and add the integrated intensity to integrated_intensity.
     *
     * This function will sum up all intensities within a window of
     * mass-to-charge. It will extract around mz +/- mz_extract_window / 2.0
     * and add the result to integrated_intensity.
     *
     * @param mz_start Start of the spectrum (m/z coordinates)
     * @param mz_it Current m/z position (will be modified)
     * @param mz_end End of the spectrum (m/z coordinates)
     * @param int_it Current intensity position (will be modified)
     * @param mz Target m/z for the current ion
     * @param integrated_intensity Resulting intensity (will be overwritten)
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param ppm Whether the parameter mz_extraction_window is given in ppm or Th
     *
     * @note This function will change the position of the iterators mz_it and
     * int_it and it can *not* extract any data if the mz-iterator is already
     * passed the mz value given. It is thus critically important to provide
     * all mz values to be extracted in ascending order!
     *
    */
    void extract_value_tophat(const std::vector<double>::const_iterator& mz_start,
                              std::vector<double>::const_iterator& mz_it,
                              const std::vector<double>::const_iterator& mz_end,
                              std::vector<double>::const_iterator& int_it,
                              const double mz,
                              double& integrated_intensity,
                              const double mz_extraction_window,
                              const bool ppm);

    /**
     * @brief Extract the next m/z value and add the integrated intensity to integrated_intensity.
     *
     * This function will sum up all intensities within a two-dimensional
     * window of mass-to-charge and ion mobility. It will extract around mz +/-
     * mz_extract_window / 2.0 and im +/- im_extraction_window / 2.0 and add
     * the result to integrated_intensity.
     *
     * @param mz_start Start of the spectrum (m/z coordinates)
     * @param mz_it Current m/z position (will be modified)
     * @param mz_end End of the spectrum (m/z coordinates)
     * @param int_it Current intensity position (will be modified)
     * @param im_it Current ion mobility position (will be modified)
     * @param mz Target m/z for the current ion
     * @param im Target ion mobility for the current ion
     * @param integrated_intensity Resulting intensity (will be overwritten)
     * @param mz_extraction_window Extracts a window of this size in m/z
     * dimension (e.g. a window of 50 ppm means an extraction of 25 ppm on
     * either side)
     * @param im_extraction_window Extracts a window of this size in ion mobility dimension.
     * @param ppm Whether the parameter mz_extraction_window is given in ppm or Th
     *
     * @note This function will change the position of the iterators mz_it,
     * int_it and im_it and it can *not* extract any data if the mz-iterator is
     * already passed the mz value given. It is thus critically important to
     * provide all mz values to be extracted in ascending order!
     *
    */
    void extract_value_tophat(const std::vector<double>::const_iterator& mz_start,
                              std::vector<double>::const_iterator& mz_it,
                              const std::vector<double>::const_iterator& mz_end,
                              std::vector<double>::const_iterator& int_it,
                              std::vector<double>::const_iterator& im_it,
                              const double mz,
                              const double im,
                              double& integrated_intensity,
                              const double mz_extraction_window,
                              const double im_extraction_window,
                              const bool ppm);

private:

    int getFilterNr_(const String& filter);

  };

}


