// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_FORMAT_DATAACCESS_SIRIUSMZTABWRITER_H
#define OPENMS_FORMAT_DATAACCESS_SIRIUSMZTABWRITER_H

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusMzTabWriter
  {
  public:

    /**
    @brief Internal structure used in @ref SiriusAdapter that is used
    for the conversion of the sirius output to an mzTab.
    @ingroup ID
    */

     /// store a specific @param number of lines from sirius output
     /// @return mzTab

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
      OpenMS::String scan_index;
      std::vector<SiriusAdapterHit> hits;
    };

    struct SiriusAdapterRun
    {
      std::vector<SiriusAdapterIdentification> identifications;
    };

    // extract scan_index from filepath
    static String extract_scan_index(const String & path);

    //Output of Sirius is one directory per spectrum/compound
    //paths: Path to output directories of sirius
    //number: Amount of entries for each file/compound should be written to the mztab file
    static void read(const std::vector<String> & paths, Size number, MzTab & result);

  };

}

#endif //OPENMS_FORMAT_DATAACCESS_SIRIUSMZTABWRITER_H
