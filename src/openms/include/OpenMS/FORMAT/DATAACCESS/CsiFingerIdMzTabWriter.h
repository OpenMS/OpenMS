// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/MzTabFile.h>

namespace OpenMS
{
  class OPENMS_DLLAPI CsiFingerIdMzTabWriter
      {
          public:

          /**
          @brief Internal structure used in @ref TOPP_SiriusAdapter that is used
           for the conversion of the Csi:FingerID output to an mzTab.

           CsiAdapterHit:
           inchikey2D (String)
           inchi (String)
           rank (int)  - Rank of the identification for a compound (spectrum) calculated by CSI:FingerID
           molecular_formula (String) - sumformula
           score (int) - Score of the identification for a compound (spectrum) calculated by CSI:FingerID
           name (String)
           smiles (String)
           pubchemids (vector<String>) - PubChemID as reference
           links (vector<String>) - Links to the database

           CsiAdapterIdentification:
           scan_index (int) - Index of the spectrum used for identification
           scan_number (int) - NativeId of the spectrum used for identification
           feature_id (String) - FeatureId (if spectrum was assigned to a feature)
           hits (vector<CsiAdapterHit>)

           CsiAdapterRun:
           identifications (vector<CSIAdapterIdentification>)

           @ingroup DATAACCESS
          */

          class CsiAdapterHit
          {
          public:
            OpenMS::String inchikey2D;
            OpenMS::String inchi;
            unsigned int rank = 0;
            unsigned int formula_rank = 0;
            OpenMS::String adduct;
            OpenMS::String molecular_formula;
            double score = 0.;
            OpenMS::String name;
            OpenMS::String smiles;
            OpenMS::String xlogp;
            OpenMS::String dbflags;
            std::vector<String> pubchemids;
            std::vector<String> links;

          };

          class CsiAdapterIdentification
          {
          public:
            double mz = 0.;
            double rt = 0.;
            OpenMS::StringList native_ids;
            int scan_index = -1;
            int scan_number = -1;
            OpenMS::String feature_id;
            std::vector<CsiAdapterHit> hits;
          };

          class CsiAdapterRun
          {
          public:
            std::vector <CsiAdapterIdentification> identifications;
          };

          /**
          @brief Conversion of CSI:FingerID output to mzTab
          
          Output of CSI:FingerID is one directory per spectrum/compound
          @param sirius_output_paths: Path to output directories of Sirius
          @param original_input_mzml: Path to original input mzml of SiriusAdapter
          @param top_n_hits: Top n  entries for each compound written to the result file
          
          @return Result written to mzTab
          */
          static void read(const std::vector<String>& sirius_output_paths,
                           const String& original_input_mzml,
                           const Size& top_n_hits,
                           MzTab& result);

      };
}
