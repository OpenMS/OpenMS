// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// forward declarations
struct sqlite3;
struct sqlite3_stmt;

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief Sqlite handler for storing spectra and chromatograms in sqMass format.

        @note Do not use this class directly, rather use SqMassFile.

        @note Due to the performance characteristics of the underlying SQLite
        database, it is highly recommended to read and write data
        (spectra/chromatograms) in batch. This is supported in this class and
        essential for reasonable performance. The current class does support
        batching SQL statements which can be controlled using setConfig and it
        is recommended to set the batch size to at least 500.
        The underlying SQLite database only stores the most essential
        parameters of a MS experiment, to store the complete meta-data, a
        zipped representation of the mzML data structure can be written
        directly into the database (and will be retrieved when converting
        back).

        This class also supports writing data using the lossy numpress
        compression format.

        This class contains the internal data structures and SQL statements for
        communication with the SQLite database

    */
    class OPENMS_DLLAPI MzMLSqliteHandler
    {

public:

      /**
          @brief Constructor of sqMass file

          @param filename The sqMass filename
          @param run_id Unique identifier which links the sqMass and OSW file. It is currently only used for storing and ignored when reading an sqMass file.
      */
      MzMLSqliteHandler(const String& filename, const UInt64 run_id);

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
      void readExperiment(MSExperiment & exp, bool meta_only = false) const;

      /// extract the RUN::ID from the sqMass file
      /// @throws Exception::SqlOperationFailed more than on run exists
      UInt64 getRunID() const;

      /**
          @brief Read an set of spectra (potentially restricted to a subset)

          @param exp The result
          @param indices A list of indices restricting the resulting spectra only to those specified here
          @param meta_only Only read the meta data
      */
      void readSpectra(std::vector<MSSpectrum> & exp, const std::vector<int> & indices, bool meta_only = false) const;

      /**
          @brief Read an set of chromatograms (potentially restricted to a subset)

          @param exp The result
          @param indices A list of indices restricting the resulting chromatograms only to those specified here
          @param meta_only Only read the meta data
      */
      void readChromatograms(std::vector<MSChromatogram> & exp, const std::vector<int> & indices, bool meta_only = false) const;

      /**
          @brief Get number of spectra in the file

          @return The number of spectra
      */
      Size getNrSpectra() const;

      /**
          @brief Get number of chromatograms in the file

          @return The number of chromatograms
      */
      Size getNrChromatograms() const;

      /**
          @brief Set file configuration

          @param write_full_meta Whether to write a complete mzML meta data structure into the RUN_EXTRA field (allows complete recovery of the input file)
          @param use_lossy_compression Whether to use lossy compression (ms numpress)
          @param linear_abs_mass_acc Accepted loss in mass accuracy (absolute m/z, in Th)
          @param sql_batch_size Batch size of SQL insert statements
      */
      void setConfig(bool write_full_meta, bool use_lossy_compression, double linear_abs_mass_acc, int sql_batch_size = 500) 
      {
        write_full_meta_ = write_full_meta;
        use_lossy_compression_ = use_lossy_compression;
        linear_abs_mass_acc_ = linear_abs_mass_acc; 
        sql_batch_size_ = sql_batch_size; 
      }

      /**
          @brief Get spectral indices around a specific retention time

          @param RT The retention time
          @param deltaRT Tolerance window around RT (if less or equal than zero, only the first spectrum *after* RT is returned)
          @param indices Spectra to consider (if empty, all spectra are considered)
          @return The indices of the spectra within RT +/- deltaRT
      */
      std::vector<size_t> getSpectraIndicesbyRT(double RT, double deltaRT, const std::vector<int> & indices) const;

protected:

      void populateChromatogramsWithData_(sqlite3 *db, std::vector<MSChromatogram>& chromatograms) const;

      void populateChromatogramsWithData_(sqlite3 *db, std::vector<MSChromatogram>& chromatograms, const std::vector<int> & indices) const;

      void populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum>& spectra) const;

      void populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum>& spectra, const std::vector<int> & indices) const;

      void prepareChroms_(sqlite3 *db, std::vector<MSChromatogram>& chromatograms, const std::vector<int> & indices = {}) const;

      void prepareSpectra_(sqlite3 *db, std::vector<MSSpectrum>& spectra, const std::vector<int> & indices = {}) const;
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

          @note Be careful with this function, calling this on an existing file
                will delete the file!
      */
      void createTables();

      /**
          @brief Writes a set of spectra to disk

          @param spectra The spectra to write
      */
      void writeSpectra(const std::vector<MSSpectrum>& spectra);

      /**
          @brief Writes a set of chromatograms to disk

          @param chroms The chromatograms to write
      */
      void writeChromatograms(const std::vector<MSChromatogram>& chroms);

      /**
          @brief Write the run-level information for an experiment into tables

          @note This is a low level function, do not call this function unless you know what you are doing!

          @param exp The result data structure
          @param write_full_meta Add full meta information into sql tables
      */
      void writeRunLevelInformation(const MSExperiment& exp, bool write_full_meta);

protected:

      void createIndices_();
      //@}

      String filename_;

      /*
       * These are spectra and chromatogram ids that are global for a specific
       * database file. Keeping track of them allows us to append spectra and
       * chromatograms multiple times to a database.
       *
       * However, currently they are initialized to zero when opening a new
       * file, so appending to an existing file won't work.
      */
      Int spec_id_;
      Int chrom_id_;
      UInt64 run_id_;

      bool use_lossy_compression_;
      double linear_abs_mass_acc_; 
      double write_full_meta_; 
      int sql_batch_size_; 
    };


  } // namespace Internal
} // namespace OpenMS


