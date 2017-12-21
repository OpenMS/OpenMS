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

#ifndef OPENMS_HOST_CSIFINGERIDMZTABWRITER_H
#define OPENMS_HOST_CSIFINGERIDMZTABWRITER_H

namespace OpenMS
{
  class OPENMS_DLLAPI CsiFingerIdMzTabWriter
      {
          public:

          /**
          @brief Internal structure used in @ref SiriusAdapter that is used
           for the conversion of the Csi:FingerID output to an mzTab.
           @ingroup DATAACCESS
          */

          struct CsiAdapterHit
          {
            OpenMS::String inchikey2D;
            OpenMS::String inchi;
            unsigned int rank;
            OpenMS::String molecular_formula;
            double score;
            OpenMS::String name;
            OpenMS::String smiles;
            std::vector<String> pubchemids;
            std::vector<String> links;

          };

          struct CsiAdapterIdentification
          {
            OpenMS::String scan_index;
            std::vector<CsiAdapterHit> hits;
          };

          struct CsiAdapterRun
          {
            std::vector <CsiAdapterIdentification> identifications;
          };

          //Output of Sirius is one directory per spectrum/compound
          //paths: Path to output directories of sirius
          //number: Amount of entries for each file/compound should be written to the mztab file
          static void read(const std::vector<String> & paths, Size number, MzTab & result);

      };
}

#endif //OPENMS_HOST_CSIFINGERIDMZTABWRITER_H
