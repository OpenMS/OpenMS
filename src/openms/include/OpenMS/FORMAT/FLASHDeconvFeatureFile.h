// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
