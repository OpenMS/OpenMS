// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <numeric>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI OPXLSpectrumProcessingAlgorithms
  {
    public:

    /**
       * @brief Merges two spectra into one while correctly considering metainfo in DataArrays
       * @param[in,out] first_spectrum
       * @param[in,out] second_spectrum
       * @return A PeakSpectrum containing all peaks from both input spectra
       */
      static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

      /**
       * @brief Preprocesses spectra
       *
       * Filters out spectra with too few peaks (based on @p peptide_min_size) and those that do not fit into the precursor charge range.
       * Removes zero intensity peaks and normalizes intensities from all spectra in @p exp.
       * MS2 spectra are reduced in peak count using a WindowMower, if @p labeled is false.
       * The number of returned spectra (MS2 only) is equal to the number of input MS2 spectra for labeled data (otherwise not necessarily).
       *
       * @param[in,out] exp Input data, which contains MS1 and MS2 spectra. Will be modified.
       * @param[in] fragment_mass_tolerance For deisotoping (used only if @p deisotope is true)
       * @param[in] fragment_mass_tolerance_unit_ppm For deisotoping (used only if @p deisotope is true)
       * @param[in] peptide_min_size Minimum length of a peptide (used to filter spectra with few peaks)
       * @param[in] min_precursor_charge MS2 spectra's minimal PC charge (ignored if @p labeled is true)
       * @param[in] max_precursor_charge MS2 spectra's maximal PC charge (ignored if @p labeled is true)
       * @param[in] deisotope Deisotope MS2 spectra?
       * @param[in] labeled Is the data labeled? (i.e. keep all MS2 spectra irrespective of precursor)
       * @return A PeakMap of preprocessed spectra (MS2 spectra only)
       */
      static PeakMap preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, bool deisotope, bool labeled);

      /**
       * @brief Computes a spectrum alignment while considering fragment charges stored in a IntegerDataArray and a cut-off for the intensity difference ratio
       * @param[out] alignment The empty alignment, that will be filled by the algorithm
       * @param[in] fragment_mass_tolerance The peak mass tolerance
       * @param[in] fragment_mass_tolerance_unit_ppm True if the given tolerance is a ppm tolerance, false if tolerance is in Da
       * @param[in] theo_spectrum The first spectrum to be aligned (preferably the theoretical one)
       * @param[in] exp_spectrum the second spectrum to be aligned (preferably the experimental one)
       * @param[in] theo_charges IntegerDataArray with charges for the theo_spectrum
       * @param[in] exp_charges IntegerDataArray with charges for the exp_spectrum
       * @param[out] ppm_error_array empty FloatDataArray to be filled with per peak ppm errors
       * @param[in] intensity_cutoff Peaks will only be aligned if intensity1 / intensity2 > intensity_cutoff, with intensity1 being the lower of the two compared peaks and intensity2 the higher one. Set to 0 to ignore intensity differences.
       */
      static void getSpectrumAlignmentFastCharge(
            std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance,
            bool fragment_mass_tolerance_unit_ppm,
            const PeakSpectrum& theo_spectrum,
            const PeakSpectrum& exp_spectrum,
            const DataArrays::IntegerDataArray& theo_charges,
            const DataArrays::IntegerDataArray& exp_charges,
            DataArrays::FloatDataArray& ppm_error_array,
            double intensity_cutoff = 0.0);

            /**
             * @brief Computes a spectrum alignment while considering fragment charges. Uses TSGXLMS::SimplePeak for the theoretical spectrum and its charges. Does not consider intensities.
             * @param[out] alignment The empty alignment, that will be filled by the algorithm
             * @param[in] fragment_mass_tolerance The peak mass tolerance
             * @param[in] fragment_mass_tolerance_unit_ppm True if the given tolerance is a ppm tolerance, false if tolerance is in Da
             * @param[in] theo_spectrum The first spectrum to be aligned (preferably the theoretical one)
             * @param[in] exp_spectrum the second spectrum to be aligned (preferably the experimental one)
             * @param[in] exp_charges IntegerDataArray with charges for the exp_spectrum
             */
            static void getSpectrumAlignmentSimple(
                  std::vector<std::pair<Size, Size> > & alignment,
                  double fragment_mass_tolerance,
                  bool fragment_mass_tolerance_unit_ppm,
                  const std::vector< SimpleTSGXLMS::SimplePeak >& theo_spectrum,
                  const PeakSpectrum& exp_spectrum,
                  const DataArrays::IntegerDataArray& exp_charges);
  };

}
