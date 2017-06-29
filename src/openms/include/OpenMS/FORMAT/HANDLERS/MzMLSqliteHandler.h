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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLSQLITEHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZMLSQLITEHANDLER_H

#include <OpenMS/KERNEL/MSExperiment.h>

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
        @brief Sqlite handler for storing spectra and chromatograms

        @note Do not use this class directly, rather use SqMassFile.

        @note This class writes spectra and chromatograms from a cache to make
        writing substantially faster. It is thus recommended to provide many
        spectra / chromatograms together to the writing function or else
        performance suffers.
    */
    class OPENMS_DLLAPI MzMLSqliteHandler
    {

public:

      MzMLSqliteHandler(String filename);

      /**@name Functions for reading files 
       *
       * ----------------------------------- 
       * Reading of SQL file starts here
       * ----------------------------------- 
       */
      //@{

      /**
          @brief Read an experiment into an MSExperiment structure

          @param exp The result data structure
          @param meta_only Only read the meta data
      */
      void readExperiment(MSExperiment & exp, bool meta_only = false);

      void readSpectra(std::vector<MSSpectrum<> > & exp, std::vector<int> & indices, bool meta_only = false);

      Size getNrSpectra() const;

      Size getNrChromatograms() const;

protected:

      void populateChromatogramsWithData_(sqlite3 *db, std::vector<MSChromatogram<> >& chromatograms);

      void populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra);

      void populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra, std::vector<int> & indices);

      void prepareChroms_(sqlite3 *db, std::vector<MSChromatogram<> >& chromatograms);

      void prepareSpectra_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra);
      //@}

public:

      /**@name Functions for writing files 
       *
       * ----------------------------------- 
       * Writing to SQL file starts here
       * ----------------------------------- 
       */
      //@{

      /**
          @brief Write an experiment to disk

          @param exp The data to write
      */
      void writeExperiment(const MSExperiment & exp);

      /**
          @brief Create data tables for a new file

          @note It is required to call this function first before writing any
                data to disk, otherwise the tables will not be set up!
      */
      void createTables();

      /**
          @brief Writes a set of spectra to disk

          @param spectra The spectra to write
      */
      void writeSpectra(const std::vector<MSSpectrum<> >& spectra);

      /**
          @brief Writes a set of chromatograms to disk

          @param chromatograms The chromatograms to write
      */
      void writeChromatograms(const std::vector<MSChromatogram<> >& chroms);

protected:

      void executeBlobBind_(sqlite3 *db, String& prepare_statement, std::vector<String>& data);

      void executeSql_(sqlite3 *db, const std::stringstream& statement);
      //@}

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

#endif // OPENMS_FORMAT_HANDLERS_MZMLSQLITEHANDLER_H

