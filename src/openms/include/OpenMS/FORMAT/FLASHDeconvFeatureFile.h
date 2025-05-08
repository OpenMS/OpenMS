// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
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
    static void writeHeader(std::fstream& fs, bool report_decoy = false);

    /// write header line for topFD feature file
    static void writeTopFDFeatureHeader(std::fstream& fs, uint ms_level);

    /// write the features in regular file output
    static void writeFeatures(const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const String& file_name, std::fstream& fs, bool report_decoy = false);

    /**
     * @brief Find mass features and write features in TopFD format files.
     * @param deconvolved_spectra deconvolved spectra - feature indices are updated only for TopFD and TopPIC outputs
     * @param mass_features mass features to be written
     * @param scan_rt_map scan number to retention time map
     * @param file_name input spectrum file name
     * @param fs file stream
     * @param ms_level ms level
     */

    static void writeTopFDFeatures(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features,
                                   const std::map<int, double>& scan_rt_map, const String& file_name, std::fstream& fs, uint ms_level);

  };
} // namespace OpenMS
