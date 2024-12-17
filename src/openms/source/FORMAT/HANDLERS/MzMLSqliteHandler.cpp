// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/MzMLFile.h> // for writing to stringstream
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <QtCore/QFileInfo>

// #include <type_traits> // for template arg detection
#include <boost/type_traits.hpp>

#include <sqlite3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>

namespace OpenMS::Internal
{

    namespace Sql = Internal::SqliteHelper;

    /*
     * @brief Helper function to concatenate integers with ","
     *
     * @param The integers to concatenate
     * 
     */
    String integerConcatenateHelper(const std::vector<int> & indices)
    {
      String tmp;
      // each element has a size of the "," character plus n digits in base10 
      tmp.reserve( int(log10(indices.size())+2) * indices.size() );
      for (Size k = 0; k < indices.size(); k++)
      {
        tmp += String(indices[k]) + ",";
      }
      tmp.resize(tmp.size() - 1); // remove last ","
      return tmp;
    }

    /*
     *
     * This function populates a set of empty data containers (MSSpectrum or
     * MSChromatogram) with data which are read from an SQLite statement. It is
     * used when reading sqMass files.  It parses all rows produced by an sql
     * statement with the following columns:
     *
     * id (integer)
     * native_id (string)
     * compression (int)
     * data_type (int)
     * binary_Data (blob)
     *
     * It is designed to work with containers of type MSSpectrum and
     * MSChromatogram to provide a single function for both use-cases.
     *
     */
    template<class ContainerT>
    void populateContainer_sub_(sqlite3_stmt *stmt, std::vector<ContainerT>& containers)
    {
      // perform first step
      sqlite3_step(stmt);

      std::vector<int> cont_data;
      cont_data.resize(containers.size());
      std::map<Size,Size> sql_container_map;
      std::vector<double> data;
      String stemp;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        Size id_orig = sqlite3_column_int( stmt, 0 );

        // map the sql table id to the index in the "containers" vector
        if (sql_container_map.find(id_orig) == sql_container_map.end())
        {
          Size tmp = sql_container_map.size();
          sql_container_map[id_orig] = tmp;
        }
        Size curr_id = sql_container_map[id_orig];

        const unsigned char * native_id_ = sqlite3_column_text(stmt, 1);
        std::string native_id(reinterpret_cast<const char*>(native_id_), sqlite3_column_bytes(stmt, 1));

        if (curr_id >= containers.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Data for non-existent spectrum / chromatogram found");
        }
        if (native_id != containers[curr_id].getNativeID())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              String("Native id for spectrum / chromatogram does not match: ") + native_id + " != " +  containers[curr_id].getNativeID() );
        }

        int compression = sqlite3_column_int( stmt, 2 );
        int data_type = sqlite3_column_int( stmt, 3 );

        const void * raw_text = sqlite3_column_blob(stmt, 4);
        size_t blob_bytes = sqlite3_column_bytes(stmt, 4);

        // data_type is one of 0 = mz, 1 = int, 2 = rt
        // compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        data.clear();
        stemp.clear();
        if (compression == 1)
        {
          OpenMS::ZlibCompression::uncompressString(raw_text, blob_bytes, stemp);

          void* byte_buffer = reinterpret_cast<void *>(&stemp[0]);
          Size buffer_size = stemp.size();
          const double* float_buffer = reinterpret_cast<const double *>(byte_buffer);
          if (buffer_size % sizeof(double) != 0)
          {
            throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount?");
          }
          Size float_count = buffer_size / sizeof(double);
          // copy values
          data.assign(float_buffer, float_buffer + float_count);
        }
        else if (compression == 5)
        {
          OpenMS::ZlibCompression::uncompressString(raw_text, blob_bytes, stemp);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("linear");
          MSNumpressCoder().decodeNPRaw(stemp, data, config);
        }
        else if (compression == 6)
        {
          OpenMS::ZlibCompression::uncompressString(raw_text, blob_bytes, stemp);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("slof");
          MSNumpressCoder().decodeNPRaw(stemp, data, config);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Compression not supported");
        }

        if (data_type == 1)
        {
          // intensity
          if (containers[curr_id].empty())
          {
            containers[curr_id].resize(data.size());
          }
          std::vector< double >::iterator data_it = data.begin();
          for (auto it = containers[curr_id].begin(); it != containers[curr_id].end(); ++it, ++data_it)
          {
            it->setIntensity(*data_it);
          }
          cont_data[curr_id] += 1;
        }
        else if (data_type == 0)
        {
          // mz (should only occur in spectra)
          if (boost::is_same<ContainerT, MSChromatogram>::value) 
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                "Found m/z data type for chromatogram (instead of retention time)");
          }

          if (containers[curr_id].empty())
          {
            containers[curr_id].resize(data.size());
          }
          std::vector< double >::iterator data_it = data.begin();
          for (auto it = containers[curr_id].begin(); it != containers[curr_id].end(); ++it, ++data_it)
          {
            it->setPos(*data_it);
          }
          cont_data[curr_id] += 1;
        }
        else if (data_type == 2)
        {
          // rt (should only occur in chromatograms)
          if (boost::is_same<ContainerT, MSSpectrum >::value) 
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                "Found retention time data type for spectrum (instead of m/z)");
          }
          if (containers[curr_id].empty()) containers[curr_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (auto it = containers[curr_id].begin(); it != containers[curr_id].end(); ++it, ++data_it)
          {
            it->setPos(*data_it);
          }
          cont_data[curr_id] += 1;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Found data type other than RT/Intensity for spectra");
        }

        sqlite3_step( stmt );
      }

      // ensure that all spectra/chromatograms have their data: we expect two data arrays per container (int and mz/rt)
      for (Size k = 0; k < cont_data.size(); k++)
      {
        if (cont_data[k] < 2)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              String("Spectrum/Chromatogram ") + k + " does not have 2 data arrays.");
        }
      }
    }

    // the cost for initialization and copy should be minimal
    //  - a single C string is created
    //  - two ints
    MzMLSqliteHandler::MzMLSqliteHandler(const String& filename, const UInt64 run_id) :
      filename_(filename),
      spec_id_(0),
      chrom_id_(0),
      run_id_(Internal::SqliteHelper::clearSignBit(run_id)),
      use_lossy_compression_(true),
      linear_abs_mass_acc_(0.0001), // set the desired mass accuracy = 1ppm at 100 m/z
      write_full_meta_(true)
    {
    }

    void MzMLSqliteHandler::readExperiment(MSExperiment& exp, bool meta_only) const
    {
      SqliteConnector conn(filename_);
      sqlite3* db = conn.getDB();

      Size nr_results = 0;
      if (write_full_meta_)
      {
        std::string select_sql = "SELECT " \
          "RUN.ID as run_id," \
          "RUN.NATIVE_ID as native_id," \
          "RUN.FILENAME as filename," \
          "RUN_EXTRA.DATA as data " \
          "FROM RUN " \
          "LEFT JOIN RUN_EXTRA ON RUN.ID = RUN_EXTRA.RUN_ID " \
          ";";

        sqlite3_stmt* stmt;
        conn.prepareStatement(&stmt, select_sql);
        sqlite3_step(stmt);

        // read data (throw exception if we find multiple runs)
        while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
        {
          if (nr_results > 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "More than one run found, cannot read both into memory");
          }

          const void* raw_text = sqlite3_column_blob(stmt, 3);
          size_t blob_bytes = sqlite3_column_bytes(stmt, 3);

          // create mzML file and parse full structure
          if (blob_bytes > 0)
          {
            MzMLFile f;
            std::string uncompressed;
            OpenMS::ZlibCompression::uncompressString(raw_text, blob_bytes, uncompressed);
            f.loadBuffer(uncompressed, exp);

            nr_results++;
          }
          else
          {
            const unsigned char* native_id = sqlite3_column_text(stmt, 1);
            const unsigned char* filename = sqlite3_column_text(stmt, 2);
            OPENMS_LOG_WARN << "Warning: no full meta data found for run " << native_id << " from file " << filename << std::endl;
          }
          sqlite3_step(stmt);
        }

        // free memory
        sqlite3_finalize(stmt);

        if (nr_results == 0)
        {
          OPENMS_LOG_WARN << "Warning: no meta data found, fall back to inference from SQL data structures." << std::endl;
        }
      }

      bool exp_empty = (exp.getNrChromatograms() == 0 && exp.getNrSpectra() == 0);
      if (!write_full_meta_ || nr_results == 0 || exp_empty)
      {
        // creates the spectra and chromatograms but does not fill them with
        // data (provides option to return meta-data only)
        std::vector<MSChromatogram> chromatograms;
        std::vector<MSSpectrum> spectra;
        prepareChroms_(db, chromatograms);
        prepareSpectra_(db, spectra);
        exp.setChromatograms(chromatograms);
        exp.setSpectra(spectra);
      }

      exp.setSqlRunID(getRunID());

      if (meta_only)
      {
        return;
      }

      populateChromatogramsWithData_(db, exp.getChromatograms());
      populateSpectraWithData_(db, exp.getSpectra());
    }

    UInt64 MzMLSqliteHandler::getRunID() const
    {
      SqliteConnector conn(filename_);
      Size nr_results = 0;
      
      std::string select_sql = "SELECT RUN.ID FROM RUN;";

      sqlite3_stmt* stmt;
      conn.prepareStatement(&stmt, select_sql);
      Sql::SqlState state = Sql::SqlState::SQL_ROW;
      UInt64 id;
      while ((state = Sql::nextRow(stmt, state)) == Sql::SqlState::SQL_ROW)
      {
        ++nr_results;
        id = Sql::extractInt64(stmt, 0);
      }
      // free memory
      sqlite3_finalize(stmt);

      if (nr_results != 1)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "File '" + filename_ + "' contains more than one run. This is currently not supported!");
      }
      return id;
    }

    void MzMLSqliteHandler::readSpectra(std::vector<MSSpectrum> & exp, const std::vector<int> & indices, bool meta_only) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index")

      // creates the spectra but does not fill them with data (provides option
      // to return meta-data only)
      SqliteConnector conn(filename_);
      prepareSpectra_(conn.getDB(), exp, indices);
      if (indices.size() != exp.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
            String("Illegal spectral indices detected ") + integerConcatenateHelper(indices) + \
            " for file of size " + getNrSpectra());
      }

      if (meta_only)
      {
        return;
      }

      populateSpectraWithData_(conn.getDB(), exp, indices);
    }

    void MzMLSqliteHandler::readChromatograms(std::vector<MSChromatogram> & exp,
                                              const std::vector<int> & indices,
                                              bool meta_only) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index")

      // creates the chromatograms but does not fill them with data (provides
      // option to return meta-data only)
      SqliteConnector conn(filename_);
      prepareChroms_(conn.getDB(), exp, indices);
      if (indices.size() != exp.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
            String("Illegal chromatogram indices detected ") + integerConcatenateHelper(indices) + \
            " for file of size " + getNrChromatograms());
      }

      if (meta_only)
      {
        return;
      }

      populateChromatogramsWithData_(conn.getDB(), exp, indices);
    }

    Size MzMLSqliteHandler::getNrSpectra() const
    {
      SqliteConnector conn(filename_);

      int ret(0);
      sqlite3_stmt* stmt;
      conn.prepareStatement(&stmt, "SELECT COUNT(*) FROM SPECTRUM;");
      sqlite3_step(stmt);
      Sql::extractValue<int>(&ret, stmt, 0);
      sqlite3_finalize(stmt);

      return (Size)ret;
    }

    std::vector<size_t> MzMLSqliteHandler::getSpectraIndicesbyRT(double RT,
                                                                 double deltaRT,
                                                                 const std::vector<int>& indices) const
    {
      // this is necessary for some applications such as the m/z correction
      SqliteConnector conn(filename_);

      String select_sql = "SELECT " \
                          "SPECTRUM.ID as spec_id " \
                          "FROM SPECTRUM ";

      if (deltaRT > 0.0)
      {
        select_sql += " WHERE RETENTION_TIME BETWEEN " + String(RT - deltaRT) + " AND " + String(RT + deltaRT);
      }
      else
      {
        select_sql += " WHERE RETENTION_TIME >= " + String(RT);
      }

      // restrict by a given set of indices
      if (!indices.empty())
      {
        select_sql += " AND SPECTRUM.ID IN (" + integerConcatenateHelper(indices) + ")";
      }

      if (deltaRT <= 0.0) {select_sql += " LIMIT 1";} // only take the first spectrum larger than RT
      select_sql += " ;";

      // Execute SQL statement
      sqlite3_stmt * stmt;
      conn.prepareStatement(&stmt, select_sql);
      sqlite3_step(stmt);

      std::vector<size_t> result;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        result.push_back( sqlite3_column_int(stmt, 0) );
        sqlite3_step(stmt);
      }
      sqlite3_finalize(stmt);

      return result;
    }

    Size MzMLSqliteHandler::getNrChromatograms() const
    {
      SqliteConnector conn(filename_);
      int ret(0);

      sqlite3_stmt* stmt;
      conn.prepareStatement(&stmt, "SELECT COUNT(*) FROM CHROMATOGRAM;");
      sqlite3_step(stmt);
      Sql::extractValue<int>(&ret, stmt, 0);
      sqlite3_finalize(stmt);

      return (Size)ret;
    }

    void MzMLSqliteHandler::populateChromatogramsWithData_(sqlite3* db, std::vector<MSChromatogram>& chromatograms) const
    {
      std::string select_sql;
      select_sql = "SELECT " \
                    "CHROMATOGRAM.ID as chrom_id," \
                    "CHROMATOGRAM.NATIVE_ID as chrom_native_id," \
                    "DATA.COMPRESSION as data_compression," \
                    "DATA.DATA_TYPE as data_type," \
                    "DATA.DATA as binary_data " \
                    "FROM CHROMATOGRAM " \
                    "INNER JOIN DATA ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID " \
                    ";";


      // Execute SQL statement
      sqlite3_stmt* stmt;
      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      populateContainer_sub_<MSChromatogram>(stmt, chromatograms);
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::populateChromatogramsWithData_(sqlite3* db,
                                                           std::vector<MSChromatogram>& chromatograms,
                                                           const std::vector<int>& indices) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index.")
      OPENMS_PRECONDITION(indices.size() == chromatograms.size(), "Chromatograms and indices need to have the same length.")

      String select_sql = "SELECT " \
                          "CHROMATOGRAM.ID as chrom_id," \
                          "CHROMATOGRAM.NATIVE_ID as chrom_native_id," \
                          "DATA.COMPRESSION as data_compression," \
                          "DATA.DATA_TYPE as data_type," \
                          "DATA.DATA as binary_data " \
                          "FROM CHROMATOGRAM " \
                          "INNER JOIN DATA ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID " \
                          "WHERE CHROMATOGRAM.ID IN (";
      select_sql += integerConcatenateHelper(indices) + ");";

      // Execute SQL statement
      sqlite3_stmt* stmt;
      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      populateContainer_sub_<MSChromatogram>(stmt, chromatograms);
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::populateSpectraWithData_(sqlite3* db, std::vector<MSSpectrum>& spectra) const
    {
      std::string select_sql;
      select_sql = "SELECT " \
                    "SPECTRUM.ID as spec_id," \
                    "SPECTRUM.NATIVE_ID as spec_native_id," \
                    "DATA.COMPRESSION as data_compression," \
                    "DATA.DATA_TYPE as data_type," \
                    "DATA.DATA as binary_data " \
                    "FROM SPECTRUM " \
                    "INNER JOIN DATA ON SPECTRUM.ID = DATA.SPECTRUM_ID " \
                    ";";

      // Execute SQL statement
      sqlite3_stmt* stmt;
      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      populateContainer_sub_<MSSpectrum>(stmt, spectra);
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::populateSpectraWithData_(sqlite3* db,
                                                     std::vector<MSSpectrum>& spectra,
                                                     const std::vector<int>& indices) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index.")
      OPENMS_PRECONDITION(indices.size() == spectra.size(), "Spectra and indices need to have the same length.")

      String select_sql = "SELECT " \
                          "SPECTRUM.ID as spec_id," \
                          "SPECTRUM.NATIVE_ID as spec_native_id," \
                          "DATA.COMPRESSION as data_compression," \
                          "DATA.DATA_TYPE as data_type," \
                          "DATA.DATA as binary_data " \
                          "FROM SPECTRUM " \
                          "INNER JOIN DATA ON SPECTRUM.ID = DATA.SPECTRUM_ID " \
                          "WHERE SPECTRUM.ID IN (";
      select_sql += integerConcatenateHelper(indices) + ");";

      // Execute SQL statement
      sqlite3_stmt* stmt;
      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      populateContainer_sub_<MSSpectrum>(stmt, spectra);
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::prepareChroms_(sqlite3* db,
                                           std::vector<MSChromatogram>& chromatograms,
                                           const std::vector<int>& indices) const
    {
      sqlite3_stmt* stmt;
      std::string select_sql = "SELECT " \
                    "CHROMATOGRAM.ID as chrom_id," \
                    "CHROMATOGRAM.NATIVE_ID as chrom_native_id," \
                    "PRECURSOR.CHARGE as precursor_charge," \
                    "PRECURSOR.DRIFT_TIME as precursor_dt," \
                    "PRECURSOR.ISOLATION_TARGET as precursor_mz," \
                    "PRECURSOR.ISOLATION_LOWER as precursor_mz_lower," \
                    "PRECURSOR.ISOLATION_UPPER as precursor_mz_upper," \
                    "PRECURSOR.PEPTIDE_SEQUENCE as precursor_seq," \
                    "PRODUCT.CHARGE as product_charge," \
                    "PRODUCT.ISOLATION_TARGET as product_mz," \
                    "PRODUCT.ISOLATION_LOWER as product_mz_lower," \
                    "PRODUCT.ISOLATION_UPPER as product_mz_upper, " \
                    "PRECURSOR.ACTIVATION_METHOD as prec_activation, " \
                    "PRECURSOR.ACTIVATION_ENERGY as prec_activation_en " \
                    "FROM CHROMATOGRAM " \
                    "INNER JOIN PRECURSOR ON CHROMATOGRAM.ID = PRECURSOR.CHROMATOGRAM_ID " \
                    "INNER JOIN PRODUCT ON CHROMATOGRAM.ID = PRODUCT.CHROMATOGRAM_ID ";

      if (!indices.empty())
      {
        select_sql += String("WHERE CHROMATOGRAM.ID IN (") + integerConcatenateHelper(indices) + ")";
      }
      select_sql += ";";

      // See https://www.sqlite.org/c3ref/column_blob.html
      // The pointers returned are valid until a type conversion occurs as
      // described above, or until sqlite3_step() or sqlite3_reset() or
      // sqlite3_finalize() is called. The memory space used to hold strings
      // and BLOBs is freed automatically. Do not pass the pointers returned
      // from sqlite3_column_blob(), sqlite3_column_text(), etc. into
      // sqlite3_free().

      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      sqlite3_step(stmt);
      String tmp;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        chromatograms.resize(chromatograms.size()+1);
        MSChromatogram& chrom = chromatograms.back();
        OpenMS::Precursor& precursor = chrom.getPrecursor();
        OpenMS::Product& product = chrom.getProduct();

        if (Sql::extractValue(&tmp, stmt, 1))
        {
          chrom.setNativeID(tmp);
        }
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL)
        {
          precursor.setCharge(sqlite3_column_int(stmt, 2));
        }
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL)
        {
          precursor.setDriftTime(sqlite3_column_double(stmt, 3));
        }
        if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) 
        {
          precursor.setMZ(sqlite3_column_double(stmt, 4));
        }
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL) 
        {
          precursor.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 5));
        }
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL)
        {
          precursor.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 6));
        }
        if (Sql::extractValue(&tmp, stmt, 7)) precursor.setMetaValue("peptide_sequence", tmp);
        // if (sqlite3_column_type(stmt, 8) != SQLITE_NULL) product.setCharge(sqlite3_column_int(stmt, 8));
        if (sqlite3_column_type(stmt, 9) != SQLITE_NULL)
        {
          product.setMZ(sqlite3_column_double(stmt, 9));
        }
        if (sqlite3_column_type(stmt, 10) != SQLITE_NULL)
        {
          product.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 10));
        }
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL)
        {
          product.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 11));
        }
        if (sqlite3_column_type(stmt, 12) != SQLITE_NULL && sqlite3_column_int(stmt, 12) != -1
            && sqlite3_column_int(stmt, 12) < static_cast<int>(OpenMS::Precursor::SIZE_OF_ACTIVATIONMETHOD))
        {
          precursor.getActivationMethods().insert(static_cast<OpenMS::Precursor::ActivationMethod>(sqlite3_column_int(stmt, 12)));
        }
        if (sqlite3_column_type(stmt, 13) != SQLITE_NULL)
        {
          precursor.setActivationEnergy(sqlite3_column_double(stmt, 13));
        }
        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::prepareSpectra_(sqlite3 *db,
                                            std::vector<MSSpectrum>& spectra,
                                            const std::vector<int> & indices) const
    {
      sqlite3_stmt * stmt;
      std::string select_sql = "SELECT " \
                    "SPECTRUM.ID as spec_id," \
                    "SPECTRUM.NATIVE_ID as spec_native_id," \
                    "SPECTRUM.MSLEVEL as spec_mslevel," \
                    "SPECTRUM.RETENTION_TIME as spec_rt," \
                    "PRECURSOR.CHARGE as precursor_charge," \
                    "PRECURSOR.DRIFT_TIME as precursor_dt," \
                    "PRECURSOR.ISOLATION_TARGET as precursor_mz," \
                    "PRECURSOR.ISOLATION_LOWER as precursor_mz_lower," \
                    "PRECURSOR.ISOLATION_UPPER as precursor_mz_upper," \
                    "PRECURSOR.PEPTIDE_SEQUENCE as precursor_seq," \
                    "PRODUCT.CHARGE as product_charge," \
                    "PRODUCT.ISOLATION_TARGET as product_mz," \
                    "PRODUCT.ISOLATION_LOWER as product_mz_lower," \
                    "PRODUCT.ISOLATION_UPPER as product_mz_upper, " \
                    "SPECTRUM.SCAN_POLARITY as spec_polarity, " \
                    "PRECURSOR.ACTIVATION_METHOD as prec_activation, " \
                    "PRECURSOR.ACTIVATION_ENERGY as prec_activation_en " \
                    "FROM SPECTRUM " \
                    "LEFT JOIN PRECURSOR ON SPECTRUM.ID = PRECURSOR.SPECTRUM_ID " \
                    "LEFT JOIN PRODUCT ON SPECTRUM.ID = PRODUCT.SPECTRUM_ID ";

      if (!indices.empty())
      {
        select_sql += String("WHERE SPECTRUM.ID IN (") + integerConcatenateHelper(indices) + ")";
      }
      select_sql += ";";

      // See https://www.sqlite.org/c3ref/column_blob.html
      // The pointers returned are valid until a type conversion occurs as
      // described above, or until sqlite3_step() or sqlite3_reset() or
      // sqlite3_finalize() is called. The memory space used to hold strings
      // and BLOBs is freed automatically. Do not pass the pointers returned
      // from sqlite3_column_blob(), sqlite3_column_text(), etc. into
      // sqlite3_free().

      SqliteConnector::prepareStatement(db, &stmt, select_sql);
      sqlite3_step(stmt);
      OpenMS::String tmp;
      while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
      {
        spectra.resize(spectra.size() + 1);
        MSSpectrum& spec = spectra.back();
        OpenMS::Precursor precursor;
        OpenMS::Product product;
        if (Sql::extractValue(&tmp, stmt, 1))
        {
          spec.setNativeID(tmp);
        }
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL)
        {
          spec.setMSLevel(sqlite3_column_int(stmt, 2));
        }
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL)
        {
          spec.setRT(sqlite3_column_double(stmt, 3));
        }
        if (sqlite3_column_type(stmt, 4) != SQLITE_NULL)
        {
          precursor.setCharge(sqlite3_column_int(stmt, 4));
        }
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL)
        {
          precursor.setDriftTime(sqlite3_column_double(stmt, 5));
        }
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL)
        {
          precursor.setMZ(sqlite3_column_double(stmt, 6));
        }
        if (sqlite3_column_type(stmt, 7) != SQLITE_NULL)
        {
          precursor.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 7));
        }
        if (sqlite3_column_type(stmt, 8) != SQLITE_NULL)
        {
          precursor.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 8));
        }
        if (Sql::extractValue(&tmp, stmt, 9))
        {
          precursor.setMetaValue("peptide_sequence", tmp);
        }
        // if (sqlite3_column_type(stmt, 10) != SQLITE_NULL) product.setCharge(sqlite3_column_int(stmt, 10));
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL)
        {
          product.setMZ(sqlite3_column_double(stmt, 11));
        }
        if (sqlite3_column_type(stmt, 12) != SQLITE_NULL)
        {
          product.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 12));
        }
        if (sqlite3_column_type(stmt, 13) != SQLITE_NULL)
        {
          product.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 13));
        }
        if (sqlite3_column_type(stmt, 14) != SQLITE_NULL) 
        {
          int pol = sqlite3_column_int(stmt, 14);
          if (pol == 0)
          {
            spec.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
          }
          else 
          {
            spec.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
          }
        }
        if (sqlite3_column_type(stmt, 15) != SQLITE_NULL && sqlite3_column_int(stmt, 15) != -1
            && sqlite3_column_int(stmt, 15) < static_cast<int>(OpenMS::Precursor::SIZE_OF_ACTIVATIONMETHOD))
        {
          precursor.getActivationMethods().insert(static_cast<OpenMS::Precursor::ActivationMethod>(sqlite3_column_int(stmt, 15)));
        }
        if (sqlite3_column_type(stmt, 16) != SQLITE_NULL)
        {
          precursor.setActivationEnergy(sqlite3_column_double(stmt, 16));
        }
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL)
        {
          spec.getPrecursors().push_back(std::move(precursor));
        }
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL)
        {
          spec.getProducts().push_back(std::move(product));
        }
        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::writeExperiment(const MSExperiment& exp)
    {
      // write run level information
      writeRunLevelInformation(exp, write_full_meta_);

      // write data
      writeChromatograms(exp.getChromatograms());
      writeSpectra(exp.getSpectra());
    }

    void MzMLSqliteHandler::writeRunLevelInformation(const MSExperiment& exp, bool write_full_meta)
    {
      SqliteConnector conn(filename_);

      // prepare streams and set required precision (default is 6 digits)
      std::stringstream insert_run_sql;
      const std::string& native_id = exp.getLoadedFilePath(); // TODO escape stuff like ' (SQL inject)
      insert_run_sql << "INSERT INTO RUN (ID, FILENAME, NATIVE_ID) VALUES (" <<
            run_id_ << ",'" << native_id << "','" << native_id << "'); ";
      conn.executeStatement("BEGIN TRANSACTION");
      conn.executeStatement(insert_run_sql.str());
      conn.executeStatement("END TRANSACTION");

      if (write_full_meta)
      {
        MSExperiment meta;

        // copy experimental settings
        meta.reserveSpaceSpectra(exp.getNrSpectra());
        meta.reserveSpaceChromatograms(exp.getNrChromatograms());
        static_cast<ExperimentalSettings &>(meta) = exp;
        for (Size k = 0; k < exp.getNrSpectra(); k++)
        {
          MSSpectrum s = exp.getSpectra()[k];
          s.clear(false);
          meta.addSpectrum(s);
        }
        for (Size k = 0; k < exp.getNrChromatograms(); k++)
        {
          MSChromatogram c = exp.getChromatograms()[k];
          c.clear(false);
          meta.addChromatogram(c);
        }
        String prepare_statement = "INSERT INTO RUN_EXTRA (RUN_ID, DATA) VALUES ";
        prepare_statement += String("(") + run_id_ + ", ?)";
        std::vector<String> data;

        std::string output;
        MzMLFile().storeBuffer(output, meta);

        // write the full metadata into the sql file (compress with zlib before)
        std::string encoded_string;
        OpenMS::ZlibCompression::compressString(output, encoded_string);
        data.emplace_back(encoded_string);
        // data.push_back(output); // in case you need to debug on the uncompressed string ...
        conn.executeBindStatement(prepare_statement, data);
      }
    }

    void MzMLSqliteHandler::createTables()
    {
      // delete file if present
      QFile file (filename_.toQString());
      file.remove();

      SqliteConnector conn(filename_);

      // Create SQL structure
      char const *create_sql =

        // data table
        //  - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        //  - data_type is one of 0 = mz, 1 = int, 2 = rt
        //  - data contains the raw (blob) data for a single data array
        "CREATE TABLE DATA(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "COMPRESSION INT," \
        "DATA_TYPE INT," \
        "DATA BLOB NOT NULL" \
        ");" \

        // spectrum table
        "CREATE TABLE SPECTRUM(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "RUN_ID INT," \
        "MSLEVEL INT NULL," \
        "RETENTION_TIME REAL NULL," \
        "SCAN_POLARITY INT NULL," \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \

        // ms-run table
        "CREATE TABLE RUN(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "FILENAME TEXT NOT NULL, " \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \

        // ms-run extra table
        "CREATE TABLE RUN_EXTRA(" \
        "RUN_ID INT," \
        "DATA BLOB NOT NULL" \
        ");" \

        // chromatogram table
        "CREATE TABLE CHROMATOGRAM(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "RUN_ID INT," \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \

        // product table
        "CREATE TABLE PRODUCT(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "CHARGE INT NULL," \
        "ISOLATION_TARGET REAL NULL," \
        "ISOLATION_LOWER REAL NULL," \
        "ISOLATION_UPPER REAL NULL" \
        ");" \

        // precursor table
        "CREATE TABLE PRECURSOR(" \
        "SPECTRUM_ID INT," \
        "CHROMATOGRAM_ID INT," \
        "CHARGE INT NULL," \
        "PEPTIDE_SEQUENCE TEXT NULL," \
        "DRIFT_TIME REAL NULL," \
        "ACTIVATION_METHOD INT NULL," \
        "ACTIVATION_ENERGY REAL NULL," \
        "ISOLATION_TARGET REAL NULL," \
        "ISOLATION_LOWER REAL NULL," \
        "ISOLATION_UPPER REAL NULL" \
        ");";

      // Execute SQL statement
      conn.executeStatement(create_sql);
      createIndices_();
    }

    void MzMLSqliteHandler::createIndices_()
    {
      // Create SQL structure
      char const *create_sql =

        // data table
        "CREATE INDEX data_chr_idx ON DATA(CHROMATOGRAM_ID);" \
        "CREATE INDEX data_sp_idx ON DATA(SPECTRUM_ID);" \

        "CREATE INDEX spec_rt_idx ON SPECTRUM(RETENTION_TIME);" \
        "CREATE INDEX spec_mslevel_idx ON SPECTRUM(MSLEVEL);" \
        "CREATE INDEX spec_run_idx ON SPECTRUM(RUN_ID);" \

        "CREATE INDEX run_extra_idx ON RUN_EXTRA(RUN_ID);" \

        "CREATE INDEX chrom_run_idx ON CHROMATOGRAM(RUN_ID);" \

        "CREATE INDEX product_chr_idx ON DATA(CHROMATOGRAM_ID);" \
        "CREATE INDEX product_sp_idx ON DATA(SPECTRUM_ID);" \

        "CREATE INDEX precursor_chr_idx ON DATA(CHROMATOGRAM_ID);" \
        "CREATE INDEX precursor_sp_idx ON DATA(SPECTRUM_ID);";

      // Execute SQL statement
      SqliteConnector conn(filename_);
      conn.executeStatement(create_sql);
    }

    void MzMLSqliteHandler::writeSpectra(const std::vector<MSSpectrum>& spectra)
    {
      // prevent writing of empty data which would throw an SQL exception
      if (spectra.empty())
      {
        return;
      }
      SqliteConnector conn(filename_);

      // prepare streams and set required precision (default is 6 digits)
      std::stringstream insert_spectra_sql;
      std::stringstream insert_precursor_sql;
      std::stringstream insert_product_sql;

      insert_spectra_sql.precision(11);
      insert_precursor_sql.precision(11);
      insert_product_sql.precision(11);

      // Encoding options
      MSNumpressCoder::NumpressConfig npconfig_mz;
      npconfig_mz.estimate_fixed_point = true; // critical
      npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_mz.setCompression("linear");
      npconfig_mz.linear_fp_mass_acc = linear_abs_mass_acc_;
      MSNumpressCoder::NumpressConfig npconfig_int;
      npconfig_int.estimate_fixed_point = true; // critical
      npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_int.setCompression("slof");

      String prepare_statement = "INSERT INTO DATA (SPECTRUM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
      std::vector<String> data;
      int sql_it = 1;

      std::vector<String> encoded_strings_mz(spectra.size());
      std::vector<String> encoded_strings_int(spectra.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize k = 0; k < (SignedSize)spectra.size(); k++)
      {
        const MSSpectrum& spec = spectra[k];

        // encode mz data (zlib or np-linear + zlib)
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getMZ();
          }

          String uncompressed_str;
          String encoded_string;
          if (use_lossy_compression_)
          {
            MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_mz);
            OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
            encoded_strings_mz[k] = encoded_string;
          }
          else
          {
            std::string str_data = std::string((const char*) (&data_to_encode[0]), data_to_encode.size() * sizeof(double));
            OpenMS::ZlibCompression::compressString(str_data, encoded_string);
            encoded_strings_mz[k] = encoded_string;
          }
        }

        // encode intensity data (zlib or np-slof + zlib)
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          if (use_lossy_compression_)
          {
            MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_int);
            OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
            encoded_strings_int[k] = encoded_string;
          }
          else
          {
            std::string str_data = std::string((const char*) (&data_to_encode[0]), data_to_encode.size() * sizeof(double));
            OpenMS::ZlibCompression::compressString(str_data, encoded_string);
            encoded_strings_int[k] = encoded_string;
          }
        }
      }

      int nr_precursors = 0;
      int nr_products = 0;
      for (Size k = 0; k < spectra.size(); k++)
      {
        const MSSpectrum& spec = spectra[k];
        int polarity = (spec.getInstrumentSettings().getPolarity() == IonSource::POSITIVE); // 1 = positive
        insert_spectra_sql << "INSERT INTO SPECTRUM(ID, RUN_ID, NATIVE_ID, MSLEVEL, RETENTION_TIME, SCAN_POLARITY) VALUES (" <<
          spec_id_ << "," <<
          run_id_ << ",'" <<
          spec.getNativeID() << "'," <<
          spec.getMSLevel() << "," <<
          spec.getRT() << "," <<
          polarity << "); ";

        if (!spec.getPrecursors().empty())
        {
          if (spec.getPrecursors().size() > 1)
          {
            std::cout << "WARNING cannot store more than first precursor" << std::endl;
          }
          if (spec.getPrecursors()[0].getActivationMethods().size() > 1)
          {
            std::cout << "WARNING cannot store more than one activation method" << std::endl;
          }

          OpenMS::Precursor prec = spec.getPrecursors()[0];
          // see src/openms/include/OpenMS/METADATA/Precursor.h for activation modes
          int activation_method = -1;
          if (!prec.getActivationMethods().empty() )
          {
            activation_method = *prec.getActivationMethods().begin();
          }
          String pepseq;
          if (prec.metaValueExists("peptide_sequence"))
          {
            pepseq = prec.getMetaValue("peptide_sequence");
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, " <<
                "ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, " <<
                "ACTIVATION_METHOD, PEPTIDE_SEQUENCE) VALUES (" << 
              spec_id_ << "," << prec.getCharge() << "," << prec.getMZ() <<
              "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
              "," << prec.getDriftTime() << 
              "," << prec.getActivationEnergy() << 
              "," << activation_method << ",'" << pepseq << "'" << "); ";
          }
          else
          {
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, " << 
              "ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD) VALUES (" <<
              spec_id_ << "," << prec.getCharge() << "," << prec.getMZ() << 
              "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() << 
              "," << prec.getDriftTime() <<
              "," << prec.getActivationEnergy() <<
              "," << activation_method << "); ";
          }
          nr_precursors++;
        }

        if (!spec.getProducts().empty())
        {
          if (spec.getProducts().size() > 1)
          {
            std::cout << "WARNING cannot store more than first product" << std::endl;
          }
          OpenMS::Product prod = spec.getProducts()[0];
          insert_product_sql << "INSERT INTO PRODUCT (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, " << 
            "ISOLATION_LOWER, ISOLATION_UPPER) VALUES (" << 
            spec_id_ << "," << 0 << "," << prod.getMZ() << 
            "," << prod.getIsolationWindowLowerOffset() << "," << prod.getIsolationWindowUpperOffset() << "); ";
          nr_products++;
        }

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode mz data (zlib or np-linear + zlib)
        {
          data.push_back(encoded_strings_mz[k]);
          if (use_lossy_compression_)
          {
            prepare_statement += String("(") + spec_id_ + ", 0, 5, ?" + sql_it++ + " ),";
          }
          else
          {
            prepare_statement += String("(") + spec_id_ + ", 0, 1, ?" + sql_it++ + " ),";
          }
        }

        // encode intensity data (zlib or np-slof + zlib)
        {
          data.push_back(encoded_strings_int[k]);
          if (use_lossy_compression_)
          {
            prepare_statement += String("(") + spec_id_ + ", 1, 6, ?" + sql_it++ + " ),";
          }
          else
          {
            prepare_statement += String("(") + spec_id_ + ", 1, 1, ?" + sql_it++ + " ),";
          }
        }
        spec_id_++;

        if (sql_it > sql_batch_size_) // flush as sqlite can only handle so many bind_blob statements
        {
          // prevent writing of empty data which would throw an SQL exception
          if (!data.empty())
          {
            prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
            conn.executeBindStatement(prepare_statement, data);
          }

          data.clear();
          prepare_statement = "INSERT INTO DATA (SPECTRUM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
          sql_it = 1;
        }

      }

      // prevent writing of empty data which would throw an SQL exception
      if (!data.empty())
      {
        prepare_statement.resize( prepare_statement.size() -1 );
        conn.executeBindStatement(prepare_statement, data);
      }

      conn.executeStatement("BEGIN TRANSACTION");
      conn.executeStatement(insert_spectra_sql.str());
      if (nr_precursors > 0)
      {
        conn.executeStatement(insert_precursor_sql.str());
      }
      if (nr_products > 0)
      {
        conn.executeStatement(insert_product_sql.str());
      }
      conn.executeStatement("END TRANSACTION");
    }

    void MzMLSqliteHandler::writeChromatograms(const std::vector<MSChromatogram >& chroms)
    {
      // prevent writing of empty data which would throw an SQL exception
      if (chroms.empty())
      {
        return;
      }
      SqliteConnector conn(filename_);

      // prepare streams and set required precision (default is 6 digits)
      std::stringstream insert_chrom_sql;
      std::stringstream insert_precursor_sql;
      std::stringstream insert_product_sql;

      insert_chrom_sql.precision(11);
      insert_precursor_sql.precision(11);
      insert_product_sql.precision(11);

      // Encoding options
      MSNumpressCoder::NumpressConfig npconfig_mz;
      npconfig_mz.estimate_fixed_point = true; // critical
      npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_mz.setCompression("linear");
      npconfig_mz.linear_fp_mass_acc = 0.05; // set the desired RT accuracy (0.05 seconds)
      MSNumpressCoder::NumpressConfig npconfig_int;
      npconfig_int.estimate_fixed_point = true; // critical
      npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_int.setCompression("slof");

      String prepare_statement = "INSERT INTO DATA (CHROMATOGRAM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
      int sql_it = 1;

      // Perform encoding in parallel
      std::vector<String> encoded_strings_rt(chroms.size());
      std::vector<String> encoded_strings_int(chroms.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize k = 0; k < (SignedSize)chroms.size(); k++)
      {
        const MSChromatogram& chrom = chroms[k];
        // encode retention time data (zlib or np-linear + zlib)
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getRT();
          }

          String uncompressed_str;
          String encoded_string;
          if (use_lossy_compression_)
          {
            MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_mz);
            OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
            encoded_strings_rt[k] = encoded_string;
          }
          else
          {
            std::string str_data = std::string((const char*) (&data_to_encode[0]), data_to_encode.size() * sizeof(double));
            OpenMS::ZlibCompression::compressString(str_data, encoded_string);
            encoded_strings_rt[k] = encoded_string;
          }
        }

        // encode intensity data (zlib or np-slof + zlib)
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          if (use_lossy_compression_)
          {
            MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_int);
            OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
            encoded_strings_int[k] = encoded_string;
          }
          else
          {
            std::string str_data = std::string((const char*) (&data_to_encode[0]), data_to_encode.size() * sizeof(double));
            OpenMS::ZlibCompression::compressString(str_data, encoded_string);
            encoded_strings_int[k] = encoded_string;
          }
        }
      }

      std::vector<String> data;
      for (Size k = 0; k < chroms.size(); k++)
      {
        const MSChromatogram& chrom = chroms[k];
        insert_chrom_sql << "INSERT INTO CHROMATOGRAM (ID, RUN_ID, NATIVE_ID) VALUES (" << chrom_id_ << "," << run_id_ << ",'" << chrom.getNativeID() << "'); ";

        OpenMS::Precursor prec = chrom.getPrecursor();
        // see src/openms/include/OpenMS/METADATA/Precursor.h for activation modes
        int activation_method = -1;
        if (!prec.getActivationMethods().empty() )
        {
          activation_method = *prec.getActivationMethods().begin();
        }
        String pepseq;
        if (prec.metaValueExists("peptide_sequence"))
        {
          pepseq = prec.getMetaValue("peptide_sequence");
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, " <<
            "ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, " << 
            "ACTIVATION_METHOD, PEPTIDE_SEQUENCE) VALUES (" << 
            chrom_id_ << "," << prec.getCharge() << "," << prec.getMZ() << 
            "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
            "," << prec.getDriftTime() << 
            "," << prec.getActivationEnergy() << 
            "," << activation_method << ",'" << pepseq << "'" << "); ";
        }
        else
        {
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, " << 
            "ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD) VALUES (" << 
            chrom_id_ << "," << prec.getCharge() << "," << prec.getMZ() << 
            "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
            "," << prec.getDriftTime() << 
            "," << prec.getActivationEnergy() << 
            "," << activation_method << "); ";
        }

        OpenMS::Product prod = chrom.getProduct();
        insert_product_sql << "INSERT INTO PRODUCT (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, " << 
          "ISOLATION_LOWER, ISOLATION_UPPER) VALUES (" << 
          chrom_id_ << "," << 0 << "," << prod.getMZ() << 
          "," << prod.getIsolationWindowLowerOffset() << "," << prod.getIsolationWindowUpperOffset() << "); ";

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode retention time data (zlib or np-linear + zlib)
        {
          data.push_back(encoded_strings_rt[k]);
          if (use_lossy_compression_)
          {
            prepare_statement += String("(") + chrom_id_ + ", 2, 5, ?" + sql_it++ + " ),";
          }
          else
          {
            prepare_statement += String("(") + chrom_id_ + ", 2, 1, ?" + sql_it++ + " ),";
          }
        }

        // encode intensity data (zlib or np-slof + zlib)
        {
          data.push_back(encoded_strings_int[k]);
          if (use_lossy_compression_)
          {
            prepare_statement += String("(") + chrom_id_ + ", 1, 6, ?" + sql_it++ + " ),";
          }
          else
          {
            prepare_statement += String("(") + chrom_id_ + ", 1, 1, ?" + sql_it++ + " ),";
          }
        }
        chrom_id_++;

        if (sql_it > sql_batch_size_) // flush as sqlite can only handle so many bind_blob statements
        {
          // prevent writing of empty data which would throw an SQL exception
          if (!data.empty())
          {
            prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
            conn.executeBindStatement(prepare_statement, data);
          }

          data.clear();
          prepare_statement = "INSERT INTO DATA (CHROMATOGRAM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
          sql_it = 1;
        }

      }

      // prevent writing of empty data which would throw an SQL exception
      if (!data.empty())
      {
        prepare_statement.resize(prepare_statement.size() -1); // remove last ","
        conn.executeBindStatement(prepare_statement, data);
      }

      conn.executeStatement("BEGIN TRANSACTION");
      conn.executeStatement(insert_chrom_sql.str());
      conn.executeStatement(insert_precursor_sql.str());
      conn.executeStatement(insert_product_sql.str());
      conn.executeStatement("END TRANSACTION");
    }

} // namespace OpenMS  // namespace Internal

