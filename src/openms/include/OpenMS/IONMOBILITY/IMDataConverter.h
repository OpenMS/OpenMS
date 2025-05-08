// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <tuple>
#include <vector>

namespace OpenMS
{
    // forward declarations
    namespace DataArrays
    {
      class FloatDataArray;
    }
    enum class DriftTimeUnit;

    /**
      @brief This class converts PeakMaps and MSSpectra from/to different IM/FAIMS storage models

    */
    class OPENMS_DLLAPI IMDataConverter
    {
    public:
      /**
        @brief Splits a PeakMap into one PeakMap per FAIMS compensation voltage

        This only works with a PeakMap that has a FAIMS compensation voltage 
        (obtained via 'spec.getDriftTime()') associated with each spectrum.
        The spectra from the original PeakMap are moved to new PeakMaps,
        so the original PeakMap is unusable afterwards.

        @param exp The PeakMap
        @return Several maps, one for each CV
        @throws Exception::MissingInformation if @p exp is not FAIMS data
      */
      static std::vector<PeakMap> splitByFAIMSCV(PeakMap&& exp);

      
      /**
        @brief Split a (TimsTOF) ion mobility frame (i.e. a spectrum concatenated from multiple spectra with different IM values) into separate spectra
   
        The input @p im_frame must have a floatDataArray where IM values are annotated. If not, an exception is thrown.

        For the output spectra, the IM value is annotated once in `spec.getDriftTime()` (there is no metadata array which contains IM values, since they are all the same).
        
        Output spectra are sorted by m/z. Ranges of the experiment are updated.

        The reverse operation is `reshapeIMFrameToSingle()`.

        @param im_frame Concatenated spectrum representing an IM frame
        @return IM frame split into spectra (one per distinct IM value), sorted by m/z, with updated ranges

        @throws Exception::MissingInformation if @p im_frame does not have IM data in floatDataArrays
      */
      static MSExperiment reshapeIMFrameToMany(MSSpectrum im_frame);

      /**
         @brief Bins the ion mobility range into discrete bins and creates a new MSExperiment for each IM bin.
         
         The IM range (of the whole @p in) is divided into equally spaced IM-bins and the bin center is the new drift time (see `spec.getDriftTime()`).
         Usually multiple spectra from an IM frame (with close IM values) fall into the same bin. These spectra are merged using SpectraMerger's block-method.
         When merging m/z peaks of two MS spectra with SpectraMerger, parameters `mz_binning_width` and `mz_binning_width_unit` and used  internally.

         To avoid artifacts at the bin borders, each bin can be extended by `bin_extension_abs` on both sides. The actual overlap between adjacent bins is thus `2*bin_extension_abs`.

         @note All MS levels are binned. If you want to bin only a specific MS level, you need to filter the input MSExperiment before calling this function.

         @param in The PeakMap containing many 'wide' IM-frame spectra (where one spectrum contains multiple IM values).
         @param number_of_IM_bins Into how many bins should the ion mobility range be sliced?
         @param bin_extension_abs How much should each bin be extended at its borders? (in absolute IM units). The actual overlap between adjacent bins is thus `2*bin_extension_abs`.
         @param mz_binning_width The width of the m/z binning window, when merging spectra of the same IM-bin (in Da or ppm, see @p mz_binning_width_unit)
         @param mz_binning_width_unit The unit of the m/z binning window (Da or ppm)
         @return One MSExperiment per IM-bin and the corresponding binning borders

         @throws Exception::InvalidValue if any spectrum in @p in is missing an IM-float data array (see IMTypes::determineIMFormat(), or MSSpectrum::containsIMData())
         @throws Exception::InvalidValue if number_of_IM_bins == 0
         @throws Exception::InvalidValue if bin_extension_abs < 0
      */
      static std::tuple < std::vector<MSExperiment>, Math::BinContainer> splitExperimentByIonMobility(MSExperiment&& in,
                                                                                                      UInt number_of_IM_bins,
                                                                                                      double bin_extension_abs,
                                                                                                      double mz_binning_width,
                                                                                                      MZ_UNITS mz_binning_width_unit);

      /**
        @brief Collapses multiple MS spectra (each with its own drift time) from the same IM-frame into a single MSSpectrum (with an IM-float data array)

        Frames are recognized by having the same RT for subsequent spectra. The IM information is taken
        from each input spectrum's .getDriftTime().
        Multiple frames are allowed.
        If the input already contains IM-frames, they are simply copied.
  
        If a spectrum does not have drift time (spec.getDriftTime()), it is simply copied to the output and ignored during the collapsing process.

        @param in The input experiment with multiple spectra per frame
        @return result The output spectra collapsed to a single spectrum per frame

        @note This requires that spectra from the same frame have the same RT ("scan start time")

        The reverse operation is `reshapeIMFrameToMany()`.

        @throws Exception::InvalidValue if any spectrum has both a single drift time AND a IM-float data array (see IMTypes::determineIMFormat(), or MSSpectrum::containsIMData())
      */
      static MSExperiment reshapeIMFrameToSingle(const MSExperiment& in);

      /**
        @brief Convert from a Unit to a CV term and annotate is as the FDA's name. This is not very accurate (since we cannot decide if its 'raw' or 'binned' IM data),
               but it allows to reconstruct the unit from the IM float-data array which is annotated with this term.
  
        <table>
        <caption>This is the mapping</caption>
        <tr><th>Unit                           <th>CV term
        <tr><td>DriftTimeUnit::MILLISECOND     <td>MS:1002816 ! mean ion mobility array
        <tr><td>DriftTimeUnit::VSSC            <td>MS:1003008 ! raw inverse reduced ion mobility array
        </table>
    
        For any other unit  (e.g. FAIMS-Compensation voltage) we throw, since the PSI CV does not 
        (and should not?) have CV terms for other IM units in ion mobility arrays.

        @param[out] fda The FDA to be annotated as an IM array
        @param[in] unit The unit of the IM measurement

        @throws Exception::InvalidValue for unsupported units
      */
      static void setIMUnit(DataArrays::FloatDataArray& fda, const DriftTimeUnit unit);

      /**
        @brief Checks if the @p fda is an ion-mobility array and if so, returns the unit (either MILLISECOND or VSSC, or NONE)
        
        The name of the @p fda should correspond to a value set by setIMUnit(), but all CV names of child terms of 
        'MS:1002893 ! ion mobility array' are accepted.
        
        <table>
        <caption>This is the current mapping (all of which return true)</caption>
        <tr><th>CV term                                             <th>Unit
        <tr><td>MS:1002816 ! mean ion mobility array                <td>DriftTimeUnit::MILLISECOND
        <tr><td>MS:1003008 ! raw inverse reduced ion mobility array <td>DriftTimeUnit::VSSC
        <tr><td>MS:1002893 ! ion mobility array **                  <td>DriftTimeUnit::NONE
        </table>
        @p **) or a child term, which is not one of the terms used above.

        @param[in] fda Input array, which is tested for its name
        @param[out] unit If @p fda is an IM array, the @p unit will contain the IM unit (undefined otherwise)
        @return True if @p fda is an IM array, false otherwise
      */
      static bool getIMUnit(const DataArrays::FloatDataArray& fda, DriftTimeUnit& unit);
    };

} //end namespace OpenMS
