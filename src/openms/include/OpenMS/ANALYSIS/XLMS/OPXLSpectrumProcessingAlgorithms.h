// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
       * @param first_spectrum
       * @param second_spectrum
       * @return A PeakSpectrum containing all peaks from both input spectra
       */
      static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

      /**
       * @brief Preprocesses spectra
       *
       * Filters out spectra with too few peaks (based on peptide_min_size) and those that do not fit into the precursor charge range.
       * Removes zero intensity peaks and normalizes intensities.
       * If the given tolerance is low enough, deisotoping is performed. Otherwise only the 500 most intense peaks are kept, if the param labeled is false.
       * The number of returned spectra is equal to the number of input spectra for labeled data (otherwise not necessarily).
       *
       * @param exp
       * @param fragment_mass_tolerance_xlinks
       * @param fragment_mass_tolerance_unit_ppm
       * @param peptide_min_size
       * @param min_precursor_charge
       * @param max_precursor_charge
       * @param labeled
       * @return A PeakMap of preprocessed spectra
       */
      static PeakMap preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, std::vector<Size>& discarded_spectra, bool deisotope, bool labeled);

      /**
       * @brief Computes a spectrum alignment while considering fragment charges stored in a IntegerDataArray and a cut-off for the intensity difference ratio
       * @param alignment The empty alignment, that will be filled by the algorithm
       * @param fragment_mass_tolerance The peak mass tolerance
       * @param fragment_mass_tolerance_unit_ppm True if the given tolerance is a ppm tolerance, false if tolerance is in Da
       * @param theo_spectrum The first spectrum to be aligned (preferably the theoretical one)
       * @param exp_spectrum the second spectrum to be aligned (preferably the experimental one)
       * @param theo_charges IntegerDataArray with charges for the theo_spectrum
       * @param exp_charges IntegerDataArray with charges for the exp_spectrum
      * @param ppm_error_array empty FloatDataArray to be filled with per peak ppm errors
       * @param intensity_cutoff Peaks will only be aligned if intensity1 / intensity2 > intensity_cutoff, with intensity1 being the lower of the two compared peaks and intensity2 the higher one. Set to 0 to ignore intensity differences.
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
             * @param alignment The empty alignment, that will be filled by the algorithm
             * @param fragment_mass_tolerance The peak mass tolerance
             * @param fragment_mass_tolerance_unit_ppm True if the given tolerance is a ppm tolerance, false if tolerance is in Da
             * @param theo_spectrum The first spectrum to be aligned (preferably the theoretical one)
             * @param exp_spectrum the second spectrum to be aligned (preferably the experimental one)
             * @param theo_charges IntegerDataArray with charges for the theo_spectrum
             * @param exp_charges IntegerDataArray with charges for the exp_spectrum
            * @param ppm_error_array empty FloatDataArray to be filled with per peak ppm errors
             * @param intensity_cutoff Peaks will only be aligned if intensity1 / intensity2 > intensity_cutoff, with intensity1 being the lower of the two compared peaks and intensity2 the higher one. Set to 0 to ignore intensity differences.
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
