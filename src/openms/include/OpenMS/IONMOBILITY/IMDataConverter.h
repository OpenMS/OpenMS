// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Eugen Netz, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
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
        associated with each spectrum.
        The spectra from the original PeakMap are moved to new PeakMaps,
        so the original PeakMap is unusable afterwards.

        @param exp The PeakMap
        @return Several maps, split by CVs
        @throws Exception::MissingInformation if @p exp is not FAIMS data
      */
      static std::vector<PeakMap> splitByFAIMSCV(PeakMap&& exp);

      
      /**
        @brief Split a (TimsTOF) ion mobility frame (i.e. a spectrum concatenated from multiple spectra with different IM values) into separate spectra
   
        The input @p im_frame must have a floatDataArray where IM values are annotated. If not, an exception is thrown.

        To get some coarser binning, choose a smaller @p number_of_bins. The default creates a new bin (=spectrum in the output) for each distinct ion mobility value.
      
        @param im_frame Concatenated spectrum representing a frame
        @param number_of_bins In how many bins should the ion mobility frame be sliced? Default(-1) assigns all peaks with identical ion-mobility values to a separate spectrum.
        @return IM frame split into multiple bins (= 1 spectrum per bin)

        @throws Exception::MissingInformation if @p im_frame does not have IM data in floatDataArrays
      */
      static MSExperiment splitByIonMobility(MSSpectrum im_frame, UInt number_of_bins = -1);

      /**
         @brief Expands all (TimsTOF) ion mobility frames in the PeakMap (i.e. all IM spectra with an IM float data array) into separate spectra. Non-IM spectra are simply copied to the result.
 
         To get some coarser custom binning, choose a smaller @p number_of_bins. The default creates a new bin (=spectrum in the output) for each distinct ion mobility value.
         For custom bins, the IM range is divided into equally spaced bins and the bin center is the new drift time.

         @param in The PeakMap containing IM-frame spectra
         @param number_of_bins In how many bins should the ion mobility frame be sliced? Default(-1) assigns all peaks with identical ion-mobility values to a separate spectrum.
         @return All IM frames split into multiple bins (= 1 spectrum per bin)
      */
      static MSExperiment splitByIonMobility(MSExperiment&& in, UInt number_of_bins = -1);

      /**
        @brief Collapses multiple MS spectra (each with its own drift time) from the same IM-frame into a single MSSpectrum (with an IM-float data array)

        Frames are recognized by having the same RT for subsequent spectra. The IM information is taken
        from each input spectrum's .getDriftTime().
        Multiple frames are allowed.
        If the input already contains IM-frames, they are simply copied.
  
        If a spectrum does not have drift time (spec.getDriftTime()), it is simply copied to the output and ignored during the collapsing process.

        @param exp The input experiment with multiple spectra per frame
        @param result The output spectra collapsed to a single spectrum per frame

        @note This requires that spectra from the same frame have the same RT ("scan start time")

        @throws Exception::InvalidValue if any spectrum has both a single drift time AND a IM-float data array (see IMTypes::determineIMFormat)
      */
      static MSExperiment collapseFramesToSingle(const MSExperiment& in);

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
