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
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

#include <sqlite3.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <QtCore/QFileInfo>

namespace OpenMS
{
  namespace Internal
  {

    static int callback(void * /* NotUsed */, int argc, char **argv, char **azColName)
    {
      int i;
      for (i=0; i<argc; i++)
      {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }
      printf("\n");
      return(0);
    }

    // the cost for initialization and copy should be minimal
    //  - a single C string is created
    //  - two ints
    MzMLSqliteHandler::MzMLSqliteHandler(String filename) :
      filename_(filename),
      spec_id_(0),
      chrom_id_(0)
    {
    }

    sqlite3* MzMLSqliteHandler::openDB() const
    {
      sqlite3 *db;
      int rc;

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Can't open database: ") + sqlite3_errmsg(db));
      }
      return db;
    }

    void MzMLSqliteHandler::readExperiment(MSExperiment & exp, bool meta_only) const
    {
      sqlite3 *db = openDB();

      // creates the spectra and chromatograms but does not fill them with data (provides option to return meta-data only)
      std::vector<MSChromatogram<> > chromatograms;
      std::vector<MSSpectrum<> > spectra;
      prepareChroms_(db, chromatograms);
      prepareSpectra_(db, spectra);
      if (meta_only) 
      {
        exp.setChromatograms(chromatograms);
        exp.setSpectra(spectra);
        return;
      }

      populateChromatogramsWithData_(db, chromatograms);
      populateSpectraWithData_(db, spectra);

      exp.setChromatograms(chromatograms);
      exp.setSpectra(spectra);

      // free up connection
      sqlite3_close(db);
    }

    void MzMLSqliteHandler::readSpectra(std::vector<MSSpectrum<> > & exp, const std::vector<int> & indices, bool meta_only) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index")

      sqlite3 *db = openDB();

      // creates the spectra but does not fill them with data (provides option to return meta-data only)
      std::vector<MSSpectrum<> > spectra;
      prepareSpectra_(db, spectra);
      for (Size k = 0; k < indices.size(); k++)
      {
        exp.push_back(spectra[indices[k]]); // TODO make more efficient
      }

      if (meta_only) 
      {
        return;
      }

      // populateChromatogramsWithData_(db, chromatograms);
      populateSpectraWithData_(db, exp, indices);

      // free up connection
      sqlite3_close(db);
    }

    Size MzMLSqliteHandler::getNrSpectra() const
    {
      sqlite3 *db = openDB();
      sqlite3_stmt * stmt;

      Size ret(0);
      std::string select_sql;
      select_sql = "SELECT COUNT(*) FROM SPECTRUM;";
      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step(stmt);
      if (sqlite3_column_type(stmt, 0) != SQLITE_NULL) ret = sqlite3_column_int(stmt, 0);

      // free memory and free up connection
      sqlite3_finalize(stmt);
      sqlite3_close(db);

      return ret;
    }

    Size MzMLSqliteHandler::getNrChromatograms() const
    {
      sqlite3 *db = openDB();
      sqlite3_stmt * stmt;

      Size ret(0);
      std::string select_sql;
      select_sql = "SELECT COUNT(*) FROM CHROMATOGRAM;";
      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );
      if (sqlite3_column_type(stmt, 0) != SQLITE_NULL) ret = sqlite3_column_int(stmt, 0);

      // free memory and free up connection
      sqlite3_finalize(stmt);
      sqlite3_close(db);

      return ret;
    }

    void MzMLSqliteHandler::populateChromatogramsWithData_(sqlite3 *db, std::vector<MSChromatogram<> >& chromatograms) const
    {
      int rc;
      sqlite3_stmt * stmt;
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
      rc = sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      if (rc != SQLITE_OK)
      {
        std::cerr << "SQL error after sqlite3_prepare" << std::endl;
        std::cerr << "Prepared statement " << select_sql << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }
      sqlite3_step(stmt);

      std::vector<int> chromdata; chromdata.resize(chromatograms.size());
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        Size chrom_id = sqlite3_column_int( stmt, 0 );
        const unsigned char * native_id_ = sqlite3_column_text(stmt, 1);
        std::string native_id(reinterpret_cast<const char*>(native_id_), sqlite3_column_bytes(stmt, 1));

        if (chrom_id >= chromatograms.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Data for non-existent chromatogram found");
        }
        if (native_id != chromatograms[chrom_id].getNativeID())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Native id for chromatogram doesnt match");
        }

        int compression = sqlite3_column_int( stmt, 2 );
        int data_type = sqlite3_column_int( stmt, 3 );

        const void * tt = sqlite3_column_blob(stmt, 4);
        size_t blob_bytes = sqlite3_column_bytes(stmt, 4);

        // data_type is one of 0 = mz, 1 = int, 2 = rt
        // compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        std::vector<double> data;
        if (compression == 5)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("linear");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else if (compression == 6)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("slof");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Compression not supported");
        }

        if (data_type == 1)
        {
          // intensity
          if (chromatograms[chrom_id].empty()) chromatograms[chrom_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSChromatogram<>::iterator it = chromatograms[chrom_id].begin(); it != chromatograms[chrom_id].end(); ++it, ++data_it)
          {
            it->setIntensity(*data_it);
          }
          chromdata[chrom_id] += 1;
        }
        else if (data_type == 2)
        {
          // rt
          if (chromatograms[chrom_id].empty()) chromatograms[chrom_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSChromatogram<>::iterator it = chromatograms[chrom_id].begin(); it != chromatograms[chrom_id].end(); ++it, ++data_it)
          {
            it->setRT(*data_it);
          }
          chromdata[chrom_id] += 1;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Found data type other than RT/Intensity for chromatograms");
        }

        sqlite3_step( stmt );
      }

      // ensure that all spectra have their data: we expect two data arrays per spectrum (int and rt)
      for (Size k = 0; k < chromdata.size(); k++)
      {
        if (chromdata[k] < 2)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Chromatogram ") + k + " does not have 2 data arrays.");
        }
      }

      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra) const
    {
      sqlite3_stmt * stmt;
      int rc;
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
      rc = sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      if (rc != SQLITE_OK)
      {
        std::cerr << "SQL error after sqlite3_prepare" << std::endl;
        std::cerr << "Prepared statement " << select_sql << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }
      sqlite3_step(stmt);

      std::vector<int> specdata; specdata.resize(spectra.size());
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        Size spec_id = sqlite3_column_int( stmt, 0 );
        const unsigned char * native_id_ = sqlite3_column_text(stmt, 1);
        std::string native_id(reinterpret_cast<const char*>(native_id_), sqlite3_column_bytes(stmt, 1));

        if (spec_id >= spectra.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Data for non-existent spectrum found");
        }
        if (native_id != spectra[spec_id].getNativeID())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Native id for spectrum doesnt match");
        }

        int compression = sqlite3_column_int( stmt, 2 );
        int data_type = sqlite3_column_int( stmt, 3 );

        const void * tt = sqlite3_column_blob(stmt, 4);
        size_t blob_bytes = sqlite3_column_bytes(stmt, 4);

        // data_type is one of 0 = mz, 1 = int, 2 = rt
        // compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        std::vector<double> data;
        if (compression == 5)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("linear");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else if (compression == 6)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("slof");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Compression not supported");
        }

        if (data_type == 1)
        {
          // intensity
          if (spectra[spec_id].empty()) spectra[spec_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSSpectrum<>::iterator it = spectra[spec_id].begin(); it != spectra[spec_id].end(); ++it, ++data_it)
          {
            it->setIntensity(*data_it);
          }
          specdata[spec_id] += 1;
        }
        else if (data_type == 0)
        {
          // mz
          if (spectra[spec_id].empty()) spectra[spec_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSSpectrum<>::iterator it = spectra[spec_id].begin(); it != spectra[spec_id].end(); ++it, ++data_it)
          {
            it->setMZ(*data_it);
          }
          specdata[spec_id] += 1;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Found data type other than RT/Intensity for spectra");
        }

        sqlite3_step( stmt );
      }

      // ensure that all spectra have their data: we expect two data arrays per spectrum (int and mz)
      for (Size k = 0; k < specdata.size(); k++)
      {
        if (specdata[k] < 2)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Spectrum ") + k + " does not have 2 data arrays.");
        }
      }

      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::populateSpectraWithData_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra, const std::vector<int> & indices) const
    {
      OPENMS_PRECONDITION(!indices.empty(), "Need to select at least one index.")
      OPENMS_PRECONDITION(indices.size() == spectra.size(), "Spectra and indices need to have the same length.")

      sqlite3_stmt * stmt;
      int rc;
      std::string select_sql;

      select_sql = "SELECT " \
                    "SPECTRUM.ID as spec_id," \
                    "SPECTRUM.NATIVE_ID as spec_native_id," \
                    "DATA.COMPRESSION as data_compression," \
                    "DATA.DATA_TYPE as data_type," \
                    "DATA.DATA as binary_data " \
                    "FROM SPECTRUM " \
                    "INNER JOIN DATA ON SPECTRUM.ID = DATA.SPECTRUM_ID " \
                    "WHERE SPECTRUM.ID IN (";

      for (Size k = 0; k < indices.size()-1; k++)
      {
        select_sql += String(indices[k]) + ",";
      }

      select_sql += String(indices[indices.size()-1]) + ");";

      std::cout << " get statement " << select_sql << std::endl;

      // Execute SQL statement
      rc = sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      if (rc != SQLITE_OK)
      {
        std::cerr << "SQL error after sqlite3_prepare" << std::endl;
        std::cerr << "Prepared statement " << select_sql << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }
      sqlite3_step(stmt);

      std::vector<int> specdata; specdata.resize(spectra.size());
      std::map<Size,Size> spec_id_map;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {

        Size spec_id_orig = sqlite3_column_int( stmt, 0 );
        // Size spec_id = 0;
        if (spec_id_map.find(spec_id_orig) == spec_id_map.end()) 
        {
          Size tmp = spec_id_map.size();
          spec_id_map[spec_id_orig] = tmp;
        }
        Size spec_id = spec_id_map[spec_id_orig];

        // std::cout << " working on spec " << spec_id_orig << " mapped to idx " << spec_id << std::endl;
        const unsigned char * native_id_ = sqlite3_column_text(stmt, 1);
        std::string native_id(reinterpret_cast<const char*>(native_id_), sqlite3_column_bytes(stmt, 1));

        if (spec_id >= spectra.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Data for non-existent spectrum found");
        }
        if (native_id != spectra[spec_id].getNativeID())
        {
          std::cout << " mismatch " << native_id << " w " <<  spectra[spec_id].getNativeID() << std::endl;
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Native id for spectrum doesnt match");
        }

        int compression = sqlite3_column_int( stmt, 2 );
        int data_type = sqlite3_column_int( stmt, 3 );

        const void * tt = sqlite3_column_blob(stmt, 4);
        size_t blob_bytes = sqlite3_column_bytes(stmt, 4);

        // data_type is one of 0 = mz, 1 = int, 2 = rt
        // compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        std::vector<double> data;
        if (compression == 5)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("linear");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else if (compression == 6)
        {
          std::string uncompressed;
          OpenMS::ZlibCompression::uncompressString(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("slof");
          MSNumpressCoder().decodeNPRaw(uncompressed, data, config);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Compression not supported");
        }

        if (data_type == 1)
        {
          // intensity
          if (spectra[spec_id].empty()) spectra[spec_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSSpectrum<>::iterator it = spectra[spec_id].begin(); it != spectra[spec_id].end(); ++it, ++data_it)
          {
            it->setIntensity(*data_it);
          }
          specdata[spec_id] += 1;
        }
        else if (data_type == 0)
        {
          // mz
          if (spectra[spec_id].empty()) spectra[spec_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSSpectrum<>::iterator it = spectra[spec_id].begin(); it != spectra[spec_id].end(); ++it, ++data_it)
          {
            it->setMZ(*data_it);
          }
          specdata[spec_id] += 1;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Found data type other than RT/Intensity for spectra");
        }

        sqlite3_step( stmt );
      }

      // ensure that all spectra have their data: we expect two data arrays per spectrum (int and mz)
      for (Size k = 0; k < specdata.size(); k++)
      {
        if (specdata[k] < 2)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Spectrum ") + k + " does not have 2 data arrays.");
        }
      }

      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::prepareChroms_(sqlite3 *db, std::vector<MSChromatogram<> >& chromatograms) const
    {
      sqlite3_stmt * stmt;
      std::string select_sql;
      select_sql = "SELECT " \
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
                    "INNER JOIN PRODUCT ON CHROMATOGRAM.ID = PRODUCT.CHROMATOGRAM_ID " \
                    ";";

      /// TODO : do we want to support reading a subset of the data (e.g. only chromatograms xx - yy)
      ///   readChromatograms_(db, stmt, chromatograms);
      /// }

      /// void readChromatograms_(sqlite3 *db, sqlite3_stmt* stmt, std::vector<MSChromatogram<> >& chromatograms)
      /// {

      // See https://www.sqlite.org/c3ref/column_blob.html
      // The pointers returned are valid until a type conversion occurs as
      // described above, or until sqlite3_step() or sqlite3_reset() or
      // sqlite3_finalize() is called. The memory space used to hold strings
      // and BLOBs is freed automatically. Do not pass the pointers returned
      // from sqlite3_column_blob(), sqlite3_column_text(), etc. into
      // sqlite3_free().

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        MSChromatogram<> chrom;

        // int chrom_id = sqlite3_column_int(stmt, 0);
        const unsigned char * native_id = sqlite3_column_text(stmt, 1);
        chrom.setNativeID( std::string(reinterpret_cast<const char*>(native_id), sqlite3_column_bytes(stmt, 1)));
        String peptide_sequence;

        OpenMS::Precursor precursor;
        OpenMS::Product product;
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) precursor.setCharge(sqlite3_column_int(stmt, 2));
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) precursor.setDriftTime(sqlite3_column_double(stmt, 3));
        if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) precursor.setMZ(sqlite3_column_double(stmt, 4));
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL) precursor.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 5));
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL) precursor.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 6));
        if (sqlite3_column_type(stmt, 7) != SQLITE_NULL) 
        {
          const unsigned char * pepseq = sqlite3_column_text(stmt, 7);
          peptide_sequence = std::string(reinterpret_cast<const char*>(pepseq), sqlite3_column_bytes(stmt, 7));
          precursor.setMetaValue("peptide_sequence", peptide_sequence);
        }
        // if (sqlite3_column_type(stmt, 8) != SQLITE_NULL) product.setCharge(sqlite3_column_int(stmt, 8));
        if (sqlite3_column_type(stmt, 9) != SQLITE_NULL) product.setMZ(sqlite3_column_double(stmt, 9));
        if (sqlite3_column_type(stmt, 10) != SQLITE_NULL) product.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 10));
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL) product.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 11));
        if (sqlite3_column_type(stmt, 12) != SQLITE_NULL && sqlite3_column_int(stmt, 12) != -1
            && sqlite3_column_int(stmt, 12) < static_cast<int>(OpenMS::Precursor::SIZE_OF_ACTIVATIONMETHOD))
        {
          precursor.getActivationMethods().insert(static_cast<OpenMS::Precursor::ActivationMethod>(sqlite3_column_int(stmt, 12)));
        }
        if (sqlite3_column_type(stmt, 13) != SQLITE_NULL) precursor.setActivationEnergy(sqlite3_column_double(stmt, 13));

        chrom.setPrecursor(precursor);
        chrom.setProduct(product);
        chromatograms.push_back(chrom);

        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);
  }

    void MzMLSqliteHandler::prepareSpectra_(sqlite3 *db, std::vector<MSSpectrum<> >& spectra) const
    {
      sqlite3_stmt * stmt;
      std::string select_sql;
      select_sql = "SELECT " \
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
                    "LEFT JOIN PRODUCT ON SPECTRUM.ID = PRODUCT.SPECTRUM_ID " \
                    ";";

      // See https://www.sqlite.org/c3ref/column_blob.html
      // The pointers returned are valid until a type conversion occurs as
      // described above, or until sqlite3_step() or sqlite3_reset() or
      // sqlite3_finalize() is called. The memory space used to hold strings
      // and BLOBs is freed automatically. Do not pass the pointers returned
      // from sqlite3_column_blob(), sqlite3_column_text(), etc. into
      // sqlite3_free().

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        MSSpectrum<> spec;

        const unsigned char * native_id = sqlite3_column_text(stmt, 1);
        spec.setNativeID( std::string(reinterpret_cast<const char*>(native_id), sqlite3_column_bytes(stmt, 1)));
        String peptide_sequence;

        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) spec.setMSLevel(sqlite3_column_int(stmt, 2));
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) spec.setRT(sqlite3_column_double(stmt, 3));

        OpenMS::Precursor precursor;
        OpenMS::Product product;
        if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) precursor.setCharge(sqlite3_column_int(stmt, 4));
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL) precursor.setDriftTime(sqlite3_column_double(stmt, 5));
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL) precursor.setMZ(sqlite3_column_double(stmt, 6));
        if (sqlite3_column_type(stmt, 7) != SQLITE_NULL) precursor.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 7));
        if (sqlite3_column_type(stmt, 8) != SQLITE_NULL) precursor.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 8));
        if (sqlite3_column_type(stmt, 9) != SQLITE_NULL) 
        {
          const unsigned char * pepseq = sqlite3_column_text(stmt, 9);
          peptide_sequence = std::string(reinterpret_cast<const char*>(pepseq), sqlite3_column_bytes(stmt, 9));
          precursor.setMetaValue("peptide_sequence", peptide_sequence);
        }
        // if (sqlite3_column_type(stmt, 10) != SQLITE_NULL) product.setCharge(sqlite3_column_int(stmt, 10));
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL) product.setMZ(sqlite3_column_double(stmt, 11));
        if (sqlite3_column_type(stmt, 12) != SQLITE_NULL) product.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 12));
        if (sqlite3_column_type(stmt, 13) != SQLITE_NULL) product.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 13));
        if (sqlite3_column_type(stmt, 14) != SQLITE_NULL) 
        {
          int pol = sqlite3_column_int(stmt, 14);
          if (pol == 0) spec.getInstrumentSettings().setPolarity(IonSource::NEGATIVE);
          else spec.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
        }
        if (sqlite3_column_type(stmt, 15) != SQLITE_NULL && sqlite3_column_int(stmt, 15) != -1
            && sqlite3_column_int(stmt, 15) < static_cast<int>(OpenMS::Precursor::SIZE_OF_ACTIVATIONMETHOD))
        {
          precursor.getActivationMethods().insert(static_cast<OpenMS::Precursor::ActivationMethod>(sqlite3_column_int(stmt, 15)));
        }
        if (sqlite3_column_type(stmt, 16) != SQLITE_NULL) precursor.setActivationEnergy(sqlite3_column_double(stmt, 16));

        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL) spec.getPrecursors().push_back(precursor);
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL) spec.getProducts().push_back(product);
        spectra.push_back(spec);

        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::writeExperiment(const MSExperiment & exp)
    {
      writeChromatograms(exp.getChromatograms());
      writeSpectra(exp.getSpectra());
    }

    void MzMLSqliteHandler::createTables()
    {
      // delete file if present
      QFile file (filename_.toQString());
      file.remove();

      sqlite3 *db = openDB();

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
        "MSLEVEL INT NULL," \
        "RETENTION_TIME REAL NULL," \
        "SCAN_POLARITY INT NULL," \
        "NATIVE_ID TEXT NOT NULL" \
        ");" \

        // chromatogram table
        "CREATE TABLE CHROMATOGRAM(" \
        "ID INT PRIMARY KEY NOT NULL," \
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
      char *zErrMsg = 0;
      int rc;
      rc = sqlite3_exec(db, create_sql, callback, 0, &zErrMsg);
      if (rc != SQLITE_OK)
      {
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
            zErrMsg);
      }
      else {
        std::cout << "Done creating tables" << std::endl;
      }
      sqlite3_close(db);
    }

    void MzMLSqliteHandler::writeSpectra(const std::vector<MSSpectrum<> >& spectra)
    {
      // prevent writing of empty data which would throw an SQL exception
      if (spectra.empty()) return;

      char *zErrMsg = 0;

      sqlite3 *db = openDB();

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
      npconfig_mz.linear_fp_mass_acc = 0.0001; // set the desired mass accuracy = 1ppm at 100 m/z
      MSNumpressCoder::NumpressConfig npconfig_int;
      npconfig_int.estimate_fixed_point = true; // critical
      npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_int.setCompression("slof");

      String prepare_statement = "INSERT INTO DATA (SPECTRUM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
      std::vector<String> data;
      int sql_it = 1;
      int nr_precursors = 0;
      int nr_products = 0;
      for (Size k = 0; k < spectra.size(); k++)
      {
        const MSSpectrum<>& spec = spectra[k];
        int polarity = (spec.getInstrumentSettings().getPolarity() == IonSource::POSITIVE); // 1 = positive
        insert_spectra_sql << "INSERT INTO SPECTRUM(ID, NATIVE_ID, MSLEVEL, RETENTION_TIME, SCAN_POLARITY) VALUES (" << 
          spec_id_ << ",'" << spec.getNativeID() << 
          "'," << spec.getMSLevel() << "," << spec.getRT() << 
          "," << polarity << "); ";

        if (!spec.getPrecursors().empty())
        {
          if (spec.getPrecursors().size() > 1) std::cout << "WARNING cannot store more than first precursor" << std::endl;
          if (spec.getPrecursors()[0].getActivationMethods().size() > 1) std::cout << "WARNING cannot store more than one activation method" << std::endl;

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
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD, PEPTIDE_SEQUENCE) VALUES (" << 
              spec_id_ << "," << prec.getCharge() << "," << prec.getMZ() <<
              "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
              "," << prec.getDriftTime() << 
              "," << prec.getActivationEnergy() << 
              "," << activation_method << ",'" << pepseq << "'" <<  "); ";
          }
          else
          {
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD) VALUES (" << 
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
          if (spec.getProducts().size() > 1) std::cout << "WARNING cannot store more than first product" << std::endl;
          OpenMS::Product prod = spec.getProducts()[0];
          insert_product_sql << "INSERT INTO PRODUCT (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER) VALUES (" << 
            spec_id_ << "," << 0 << "," << prod.getMZ() << 
            "," << prod.getIsolationWindowLowerOffset() << "," << prod.getIsolationWindowUpperOffset() << "); ";
          nr_products++;
        }

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode mz data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getMZ();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_mz);
          OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
          data.push_back(encoded_string);
          prepare_statement += String("(") + spec_id_ + ", 0, 5, ?" + sql_it++ + " ),";
        }

        // encode intensity data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_int);
          OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
          data.push_back( encoded_string );
          prepare_statement += String("(") + spec_id_ + ", 1, 6, ?" + sql_it++ + " ),";
        }
        spec_id_++;

        if (sql_it > 500) // flush after 500 spectra as sqlite can only handle so many bind_blob statments
        {
          // prevent writing of empty data which would throw an SQL exception
          if (!data.empty())
          {
            prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
            executeBlobBind_(db, prepare_statement, data);
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
        executeBlobBind_(db, prepare_statement, data);
      }

      sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);

      executeSql_(db, insert_spectra_sql);
      if (nr_precursors > 0) executeSql_(db, insert_precursor_sql);
      if (nr_products > 0) executeSql_(db, insert_product_sql);

      sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);

      sqlite3_close(db);
    }

    void MzMLSqliteHandler::writeChromatograms(const std::vector<MSChromatogram<> >& chroms)
    {
      // prevent writing of empty data which would throw an SQL exception
      if (chroms.empty()) return;

      char *zErrMsg = 0;

      sqlite3 *db = openDB();

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
      std::vector<String> data;
      int sql_it = 1;
      for (Size k = 0; k < chroms.size(); k++)
      {
        const MSChromatogram<>& chrom = chroms[k];
        insert_chrom_sql << "INSERT INTO CHROMATOGRAM (ID, NATIVE_ID) VALUES (" << chrom_id_ << ",'" << chrom.getNativeID() << "'); ";

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
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD, PEPTIDE_SEQUENCE) VALUES (" << 
            chrom_id_ << "," << prec.getCharge() << "," << prec.getMZ() << 
            "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
            "," << prec.getDriftTime() << 
            "," << prec.getActivationEnergy() << 
            "," << activation_method << ",'" << pepseq << "'" <<  "); ";
        }
        else
        {
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER, DRIFT_TIME, ACTIVATION_ENERGY, ACTIVATION_METHOD) VALUES (" << 
            chrom_id_ << "," << prec.getCharge() << "," << prec.getMZ() << 
            "," << prec.getIsolationWindowLowerOffset() << "," << prec.getIsolationWindowUpperOffset() <<
            "," << prec.getDriftTime() << 
            "," << prec.getActivationEnergy() << 
            "," << activation_method << "); ";
        }

        OpenMS::Product prod = chrom.getProduct();
        insert_product_sql << "INSERT INTO PRODUCT (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, ISOLATION_LOWER, ISOLATION_UPPER) VALUES (" << 
          chrom_id_ << "," << 0 << "," << prod.getMZ() << 
          "," << prod.getIsolationWindowLowerOffset() << "," << prod.getIsolationWindowUpperOffset() << "); ";

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode retention time data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getRT();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_mz);
          OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
          data.push_back(encoded_string);
          prepare_statement += String("(") + chrom_id_ + ", 2, 5, ?" + sql_it++ + " ),";
        }

        // encode intensity data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder().encodeNPRaw(data_to_encode, uncompressed_str, npconfig_int);
          OpenMS::ZlibCompression::compressString(uncompressed_str, encoded_string);
          data.push_back( encoded_string );
          prepare_statement += String("(") + chrom_id_ + ", 1, 6, ?" + sql_it++ + " ),";
        }
        chrom_id_++;

        if (sql_it > 500) // flush after 500 chromatograms as sqlite can only handle so many bind_blob statments
        {
          // prevent writing of empty data which would throw an SQL exception
          if (!data.empty())
          {
            prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
            executeBlobBind_(db, prepare_statement, data);
          }

          data.clear();
          prepare_statement = "INSERT INTO DATA (CHROMATOGRAM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
          sql_it = 1;
        }

      }

      // prevent writing of empty data which would throw an SQL exception
      if (!data.empty())
      {
        prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
        executeBlobBind_(db, prepare_statement, data);
      }

      sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);

      executeSql_(db, insert_chrom_sql);
      executeSql_(db, insert_precursor_sql);
      executeSql_(db, insert_product_sql);

      sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);

      sqlite3_close(db);
    }

    void MzMLSqliteHandler::executeBlobBind_(sqlite3 *db, String& prepare_statement, std::vector<String>& data)
    {
      int rc;

      // The calling procedure is responsible for deleting the compiled SQL statement using sqlite3_finalize() after it has finished with it.
      sqlite3_stmt *stmt = NULL;
      const char *curr_loc;
      rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), &stmt, &curr_loc);
      if (rc != SQLITE_OK)
      {
        std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }

      for (Size k = 0; k < data.size(); k++)
      {
        // Fifth argument is a destructor for the blob.
        // SQLITE_STATIC because the statement is finalized
        // before the buffer is freed:
        rc = sqlite3_bind_blob(stmt, k+1, data[k].c_str(), data[k].size(), SQLITE_STATIC);
        if (rc != SQLITE_OK)
        {
          std::cerr << "SQL error after sqlite3_bind_blob at iteration " << k << std::endl;
          std::cerr << "Prepared statement " << prepare_statement << std::endl;
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
        } 
      }

      rc = sqlite3_step(stmt);
      if (rc != SQLITE_DONE)
      {
        std::cerr << "SQL error after sqlite3_step" << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }

      // free memory again
      sqlite3_finalize(stmt);
    }

    void MzMLSqliteHandler::executeSql_(sqlite3 *db, const std::stringstream& statement)
    {
      char *zErrMsg = 0;
      std::string insert_str = statement.str();
      int rc = sqlite3_exec(db, insert_str.c_str(), callback, 0, &zErrMsg);
      if (rc != SQLITE_OK)
      {
        std::cerr << "Error message after sqlite3_exec" << std::endl;
        std::cerr << "Prepared statement " << statement.str() << std::endl;
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, zErrMsg);
      }
    }

  } // namespace Internal
} // namespace OpenMS

