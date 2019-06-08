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

#include <OpenMS/FORMAT/MzTabFile.h>

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusMzTabWriter
  {
  public:

    /**
    @brief Internal structure used in @ref SiriusAdapter that is used
    for the conversion of the sirius output to an mzTab.
    @ingroup ID

    SiriusAdapterHit:
    formula (String) - Sumformula
    adduct (String) - Assigned adduct
    rank (int)  - Rank of the possible sumformula for a compound (spectrum) calculated by Sirius
    score (double) - Overall score of the possible sumformula for a compound (spectrum) calculated by Sirius
    treescore (double) - Fragmentation pattern score
    isoscore (double) - Isotope pattern score
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

    struct SiriusAdapterHit
    {
      OpenMS::String formula;
      OpenMS::String adduct;
      int rank;
      double score;
      double treescore;
      double isoscore;
      int explainedpeaks;
      double explainedintensity;
    };

    struct SiriusAdapterIdentification
    {
      double mz;
      double rt;
      OpenMS::String native_id;
      int scan_index;
      int scan_number;
      OpenMS::String feature_id;
      std::vector<SiriusAdapterHit> hits;
    };

    struct SiriusAdapterRun
    {
      std::vector<SiriusAdapterIdentification> identifications;
    };

    /**
    @brief Extract scan_index from filepath
    */
    static int extract_scan_index(const String & path);
 
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

}

