// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/CsvFile.h>

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusMzTabWriter
  {
  public:

    /**
    @brief Internal structure used in @ref TOPP_SiriusAdapter that is used
    for the conversion of the sirius output to an mzTab.
    @ingroup ID

    SiriusAdapterHit:
    formula (String) - Sumformula
    adduct (String) - Assigned adduct
    precursor_formula (String) - Sumformula of the precursor (can be the same as formula)
    rank (int)  - Rank of the possible sumformula for a compound (spectrum) calculated by Sirius
    rank_score (double) - Ranking score
    iso_score (double) - Isotope pattern score
    tree_score (double) - Fragmentation pattern score
    sirius_score (double) - Overall score of the possible sumformula for a compound (spectrum) calculated by Sirius
    explainedpeaks (int) - Number of explained peaks
    explainedintensity (double) - Relative amount of explained intensity

    SiriusAdapterIdentification:
    scan_index (int) - Index of the spectrum used
    scan_number (int) - NativeId of the spectrum used
    feature_id (String) - FeatureId (if spectrum was assigned to a feature)
    hits (vector<SiriusAdapterHit>)

    SiriusAdapterRun:
    identifications (vector<SiriusAdapterIdentification>)

    Store a specific @param number of lines from sirius output
    @return mzTab
    */

    class SiriusAdapterHit
    {
    public:
      OpenMS::String formula;
      OpenMS::String adduct;
      OpenMS::String precursor_formula;
      int rank = 0;
      double iso_score = 0.0;
      double tree_score = 0.0;
      double sirius_score = 0.0;
      int explainedpeaks = 0;
      double explainedintensity = 0.0;
      double median_mass_error_fragment_peaks_ppm = 0.0;
      double median_absolute_mass_error_fragment_peaks_ppm = 0.0;
      double mass_error_precursor_ppm = 0.0;
    };

    class SiriusAdapterIdentification
    {
    public:
      double mz = 0.;
      double rt = 0.;
      OpenMS::StringList native_ids;
      int scan_index = -1;
      int scan_number = -1;
      OpenMS::String feature_id;
      std::vector<SiriusAdapterHit> hits;
    };

    class SiriusAdapterRun
    {
    public:
      std::vector<SiriusAdapterIdentification> identifications;
    };

    class SiriusSpectrumMSInfo
    {
    public:
      StringList ext_n_id; // multiple possible MS2 spectra
      double ext_mz = 0.0;
      double ext_rt = 0.0;
    };

    /**
    @brief Extract scan_index from filepath
    */
    static int extractScanIndex(const String& path);

    /**
    @brief Extract scan_number from filepath
    */
    static int extractScanNumber(const String& path);

    /**
    @brief Extract feature_id from filepath
    */
    static String extractFeatureId(const String& path);

    /**
    @brief Extract columnname and index based in SIRIUS entries
    */
    static std::map< std::string, Size > extract_columnname_to_columnindex(CsvFile& csvfile);

    /**
     @brief Extract mz, rt of the precursor and the nativeID of the corresponding MS2 spectra in the spectrum.ms file
    */
    static SiriusSpectrumMSInfo extractSpectrumMSInfo(const String& single_sirius_path);

    /**
    @brief Conversion of sirius output to mzTab
    
    Output of Sirius is one directory per spectrum/compound
    @param sirius_output_paths: Path to output directories of Sirius
    @param original_input_mzml: Path to mzml input of SiriusAdapter
    @param top_n_hits: Top n  entries for each compound written to the result file     
    
    @return: Result written to mzTab
    */
    static void read(const std::vector<String>& sirius_output_paths,
                     const String& original_input_mzml,
                     const Size& top_n_hits,
                     MzTab& result);

  };

  namespace SiriusVersion
  {
    /// SIRIUS version
    inline const std::string CURRENT_VERSION = "5.6.3";
  }

}

