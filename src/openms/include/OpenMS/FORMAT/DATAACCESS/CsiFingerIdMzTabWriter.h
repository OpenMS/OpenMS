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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  class OPENMS_DLLAPI CsiFingerIdMzTabWriter
      {
          public:

          /**
          @brief Internal structure used in @ref SiriusAdapter that is used
           for the conversion of the Csi:FingerID output to an mzTab.

           CsiAdapterHit:
           inchikey2D (String)
           inchi (String)
           rank (int)  - Rank of the identification for a compound (spectrum) calculated by CSI:FingerID
           molecular_formula (String) - sumformula
           score (int) - Score of the identification for a compound (spectrum) calculated by CSI:FingerID
           name (String)
           smiles (String)
           pubchemids (vector<String>) - Pubchemid as reference
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

          struct CsiAdapterHit
          {
            OpenMS::String inchikey2D;
            OpenMS::String inchi;
            unsigned int rank = 0;
            OpenMS::String molecular_formula;
            double score = 0.;
            OpenMS::String name;
            OpenMS::String smiles;
            std::vector<String> pubchemids;
            std::vector<String> links;

          };

          struct CsiAdapterIdentification
          {
            double mz = 0.;
            double rt = 0.;
            OpenMS::String native_id;
            int scan_index = -1;
            int scan_number = -1;
            OpenMS::String feature_id;
            std::vector<CsiAdapterHit> hits;
          };

          struct CsiAdapterRun
          {
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
          static void read(const std::vector<String> & sirius_output_paths,
                           const String & original_input_mzml,
                           const Size & top_n_hits,
                           MzTab & result);

      };
}

