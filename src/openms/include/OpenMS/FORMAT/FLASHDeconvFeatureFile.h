// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/config.h>
#include <iomanip>
#include <iostream>


namespace OpenMS
{
  /**
    @brief FLASHDeconv feature level output *.tsv, *.ms1ft (for Promex), *.feature (for TopPIC) file formats
     @ingroup FileIO
**/

  class OPENMS_DLLAPI FLASHDeconvFeatureFile
  {
  public:
    /// write header line for regular file output
    static void writeHeader(std::fstream& fs);

    /// write header line for promex file output
    static void writePromexHeader(std::fstream& fs);

    /// write header line for topFD feature files
    static void writeTopFDFeatureHeader(std::vector<std::fstream>& fs);

    /// write the features in regular file output
    static void writeFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const String& file_name, std::fstream& fs);


    /**
     * @brief Find mass features and write features in TopFD format files.
     * @param mass_features mass features to be written
     * @param precursor_peak_groups precursor peak groups of MSn spectra that are used only when topfd_feature_out is set
     * @param scan_rt_map scan number to retention time map
     * @param file_name input spectrum file name
     * @param fs file stream
     */

    static void writeTopFDFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const std::map<int, PeakGroup>& precursor_peak_groups,
                                   const std::map<int, double>& scan_rt_map, const String& file_name, std::vector<std::fstream>& fs);

    /**
     * @brief Find mass features and write features in Promex format files.
     * @param mass_features mass features to be written
     * @param precursor_peak_groups precursor peak groups of MSn spectra that are used only when topfd_feature_out is set
     * @param scan_rt_map scan number to retention time map
     * @param avg averagine to determine isotope pattern range
     * @param fs file stream
     */

    static void writePromexFeatures(const std::vector<FLASHDeconvHelperStructs::MassFeature>& mass_features, const std::map<int, PeakGroup>& precursor_peak_groups,
                                    const std::map<int, double>& scan_rt_map, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, std::fstream& fs);
  };
} // namespace OpenMS
