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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZMLSQLITESWATHHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZMLSQLITESWATHHANDLER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// forward declarations
struct sqlite3;

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief Sqlite handler for SWATH data sets

        This class represents a single sqMass file acquired in SWATH / DIA mode
        and provides some useful access to the indices of the individual SWATH
        windows.

    */
    class OPENMS_DLLAPI MzMLSqliteSwathHandler
    {

public:

      /**
          @brief Constructor

          @param filename The sqMass filename
      */
      MzMLSqliteSwathHandler(String filename) :
        filename_(filename)
      {}

      /**
          @brief Read SWATH windows boundaries from file

          @return A vector populated with SwathMap, with the following attributes initialized: center, lower and upper

      */
      std::vector<OpenSwath::SwathMap> readSwathWindows();

      /**
          @brief Read indices of MS1 spectra from file

          @return A list of spectral indices for the MS1 spectra
      */
      std::vector<int> readMS1Spectra();

      /**
          @brief Read indices of spectra belonging to specified SWATH window from file

          @param swath_map Contains the upper/lower boundaries of the SWATH window
          @return A list of spectral indices for the provided SWATH window

      */
      std::vector<int> readSpectraForWindow(OpenSwath::SwathMap swath_map);

protected:

      sqlite3* openDB();

      String filename_;

      /// Decoder/Encoder for Base64-data in MzML
      Base64 base64coder_;
      MSNumpressCoder numpress_coder_;

      /*
       * These are spectra and chromatogram ids that are global for a specific
       * database file. Keeping track of them allows us to append spectra and
       * chromatograms multiple times to a database.
      */
      Int spec_id_;
      Int chrom_id_;

    };


  }   // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_MZMLSQLITESWATHHANDLER_H
