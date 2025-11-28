// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <zlib.h>

#include <fstream>
#include <vector>
#include <cstring>

namespace OpenMS
{

  // ============================================================================
  // ZIP Archive Utilities for mzPeak format
  // ============================================================================
  //
  // The mzPeak format uses uncompressed ZIP archives (store method only) to
  // bundle multiple Parquet files. We implement minimal ZIP handling here
  // because:
  //
  // 1. mzPeak requires uncompressed (store method) ZIP only - no deflate needed
  // 2. The ZIP format for store method is straightforward binary structures
  // 3. Avoids adding libzip/minizip as an external dependency
  // 4. zlib (already required) provides CRC32 but not ZIP archive handling
  //
  // ============================================================================
  
  namespace Internal
  {
    namespace ZipFormat
    {
      // ZIP format magic signatures
      constexpr uint32_t LOCAL_FILE_HEADER_SIG = 0x04034b50;
      constexpr uint32_t CENTRAL_DIR_HEADER_SIG = 0x02014b50;
      constexpr uint32_t END_OF_CENTRAL_DIR_SIG = 0x06054b50;
      
      // Compression methods
      constexpr uint16_t COMPRESSION_STORE = 0;  // No compression (required by mzPeak)

      #pragma pack(push, 1)
      
      /// ZIP local file header (appears before each file's data)
      struct LocalFileHeader
      {
        uint32_t signature;         // LOCAL_FILE_HEADER_SIG
        uint16_t version_needed;    // Version needed to extract
        uint16_t flags;             // General purpose bit flag
        uint16_t compression;       // Compression method (0 = store)
        uint16_t mod_time;          // Last mod file time
        uint16_t mod_date;          // Last mod file date
        uint32_t crc32;             // CRC-32 checksum
        uint32_t compressed_size;   // Compressed size
        uint32_t uncompressed_size; // Uncompressed size
        uint16_t filename_length;   // File name length
        uint16_t extra_length;      // Extra field length
        // Followed by: filename (variable), extra field (variable), file data
      };

      /// ZIP central directory file header
      struct CentralDirHeader
      {
        uint32_t signature;           // CENTRAL_DIR_HEADER_SIG
        uint16_t version_made;        // Version made by
        uint16_t version_needed;      // Version needed to extract
        uint16_t flags;               // General purpose bit flag
        uint16_t compression;         // Compression method
        uint16_t mod_time;            // Last mod file time
        uint16_t mod_date;            // Last mod file date
        uint32_t crc32;               // CRC-32 checksum
        uint32_t compressed_size;     // Compressed size
        uint32_t uncompressed_size;   // Uncompressed size
        uint16_t filename_length;     // File name length
        uint16_t extra_length;        // Extra field length
        uint16_t comment_length;      // File comment length
        uint16_t disk_start;          // Disk number start
        uint16_t internal_attr;       // Internal file attributes
        uint32_t external_attr;       // External file attributes
        uint32_t local_header_offset; // Offset to local header
        // Followed by: filename (variable), extra field (variable), comment (variable)
      };

      /// ZIP end of central directory record
      struct EndOfCentralDir
      {
        uint32_t signature;           // END_OF_CENTRAL_DIR_SIG
        uint16_t disk_number;         // Number of this disk
        uint16_t disk_start;          // Disk where central directory starts
        uint16_t num_entries_disk;    // Number of entries on this disk
        uint16_t num_entries_total;   // Total number of entries
        uint32_t central_dir_size;    // Size of central directory
        uint32_t central_dir_offset;  // Offset to central directory
        uint16_t comment_length;      // ZIP file comment length
        // Followed by: comment (variable)
      };
      
      #pragma pack(pop)

    } // namespace ZipFormat

    /// Entry metadata for a file within a ZIP archive
    struct ZipEntry
    {
      std::string filename;
      uint32_t local_header_offset;
      uint32_t compressed_size;
      uint32_t uncompressed_size;
    };

    /**
      @brief Read the central directory index from a ZIP archive
      
      @param archive_path Path to the ZIP file
      @return Vector of ZipEntry describing files in the archive
      @throws Exception::FileNotFound if file cannot be opened
      @throws Exception::ParseError if ZIP structure is invalid
    */
    std::vector<ZipEntry> readZipIndex(const std::string& archive_path)
    {
      std::vector<ZipEntry> entries;
      std::ifstream file(archive_path, std::ios::binary);
      
      if (!file.is_open())
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path);
      }

      // Get file size
      file.seekg(0, std::ios::end);
      const auto file_size = file.tellg();
      
      // Search backwards for End of Central Directory signature
      // EOCD is at least 22 bytes, and comment can be up to 65535 bytes
      const size_t max_search = std::min(static_cast<size_t>(file_size), size_t(65536 + 22));
      
      file.seekg(-static_cast<std::streamoff>(max_search), std::ios::end);
      std::vector<char> buffer(max_search);
      file.read(buffer.data(), max_search);
      
      // Find EOCD signature (search from end)
      int64_t eocd_pos = -1;
      for (size_t i = max_search - 22; ; --i)
      {
        uint32_t sig;
        std::memcpy(&sig, &buffer[i], sizeof(sig));
        if (sig == ZipFormat::END_OF_CENTRAL_DIR_SIG)
        {
          eocd_pos = static_cast<int64_t>(file_size) - static_cast<int64_t>(max_search) + static_cast<int64_t>(i);
          break;
        }
        if (i == 0) break;
      }
      
      if (eocd_pos < 0)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          archive_path, "Not a valid ZIP file: End of Central Directory not found");
      }

      // Read EOCD
      file.seekg(eocd_pos);
      ZipFormat::EndOfCentralDir eocd;
      file.read(reinterpret_cast<char*>(&eocd), sizeof(eocd));
      
      // Read central directory entries
      file.seekg(eocd.central_dir_offset);
      entries.reserve(eocd.num_entries_total);
      
      for (uint16_t i = 0; i < eocd.num_entries_total; ++i)
      {
        ZipFormat::CentralDirHeader header;
        file.read(reinterpret_cast<char*>(&header), sizeof(header));
        
        if (header.signature != ZipFormat::CENTRAL_DIR_HEADER_SIG)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            archive_path, "Invalid central directory entry at index " + std::to_string(i));
        }
        
        // Read filename
        std::string filename(header.filename_length, '\0');
        file.read(&filename[0], header.filename_length);
        
        // Skip extra field and comment
        file.seekg(header.extra_length + header.comment_length, std::ios::cur);
        
        entries.push_back({
          std::move(filename),
          header.local_header_offset,
          header.compressed_size,
          header.uncompressed_size
        });
      }
      
      return entries;
    }

    /**
      @brief Extract a file from a ZIP archive into memory
      
      @param archive_path Path to the ZIP file
      @param entry ZipEntry describing the file to extract
      @return Buffer containing the extracted file data
      @throws Exception::FileNotFound if archive cannot be opened
      @throws Exception::ParseError if entry is compressed (mzPeak requires store method)
    */
    std::vector<uint8_t> extractZipEntry(const std::string& archive_path, const ZipEntry& entry)
    {
      std::ifstream file(archive_path, std::ios::binary);
      
      if (!file.is_open())
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path);
      }

      // Read local file header
      file.seekg(entry.local_header_offset);
      ZipFormat::LocalFileHeader header;
      file.read(reinterpret_cast<char*>(&header), sizeof(header));
      
      if (header.signature != ZipFormat::LOCAL_FILE_HEADER_SIG)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          archive_path, "Invalid local file header for: " + entry.filename);
      }
      
      // mzPeak format requires uncompressed (store method) ZIP
      if (header.compression != ZipFormat::COMPRESSION_STORE)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          archive_path, "mzPeak requires uncompressed ZIP entries. Entry '" + 
          entry.filename + "' uses compression method " + std::to_string(header.compression));
      }
      
      // Skip to data (past filename and extra field in local header)
      file.seekg(header.filename_length + header.extra_length, std::ios::cur);
      
      // Read file data
      std::vector<uint8_t> data(entry.uncompressed_size);
      file.read(reinterpret_cast<char*>(data.data()), entry.uncompressed_size);
      
      return data;
    }

    /**
      @brief Writer for creating uncompressed ZIP archives
      
      Creates ZIP files using the store method (no compression) as required
      by the mzPeak format specification. Uses zlib for CRC32 calculation.
    */
    class ZipWriter
    {
    public:
      explicit ZipWriter(const std::string& archive_path) 
        : file_(archive_path, std::ios::binary | std::ios::trunc)
      {
        if (!file_.is_open())
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path);
        }
      }
      
      ~ZipWriter()
      {
        if (file_.is_open())
        {
          finalize();
        }
      }
      
      /// Add a file entry to the archive
      void addEntry(const std::string& filename, const std::vector<uint8_t>& data)
      {
        // Calculate CRC32 using zlib
        uLong crc = crc32(0L, Z_NULL, 0);
        crc = crc32(crc, data.data(), static_cast<uInt>(data.size()));
        
        // Record for central directory
        entries_.push_back({
          filename,
          static_cast<uint32_t>(file_.tellp()),
          static_cast<uint32_t>(data.size()),
          static_cast<uint32_t>(crc)
        });
        
        // Write local file header
        ZipFormat::LocalFileHeader header{};
        header.signature = ZipFormat::LOCAL_FILE_HEADER_SIG;
        header.version_needed = 20;  // 2.0
        header.compression = ZipFormat::COMPRESSION_STORE;
        header.crc32 = static_cast<uint32_t>(crc);
        header.compressed_size = static_cast<uint32_t>(data.size());
        header.uncompressed_size = static_cast<uint32_t>(data.size());
        header.filename_length = static_cast<uint16_t>(filename.size());
        
        file_.write(reinterpret_cast<const char*>(&header), sizeof(header));
        file_.write(filename.data(), filename.size());
        file_.write(reinterpret_cast<const char*>(data.data()), data.size());
      }
      
      /// Finalize the archive (write central directory)
      void finalize()
      {
        if (!file_.is_open()) return;
        
        const uint32_t central_dir_offset = static_cast<uint32_t>(file_.tellp());
        
        // Write central directory headers
        for (const auto& entry : entries_)
        {
          ZipFormat::CentralDirHeader header{};
          header.signature = ZipFormat::CENTRAL_DIR_HEADER_SIG;
          header.version_made = 20;
          header.version_needed = 20;
          header.compression = ZipFormat::COMPRESSION_STORE;
          header.crc32 = entry.crc32;
          header.compressed_size = entry.size;
          header.uncompressed_size = entry.size;
          header.filename_length = static_cast<uint16_t>(entry.filename.size());
          header.local_header_offset = entry.offset;
          
          file_.write(reinterpret_cast<const char*>(&header), sizeof(header));
          file_.write(entry.filename.data(), entry.filename.size());
        }
        
        const uint32_t central_dir_size = static_cast<uint32_t>(file_.tellp()) - central_dir_offset;
        
        // Write end of central directory
        ZipFormat::EndOfCentralDir eocd{};
        eocd.signature = ZipFormat::END_OF_CENTRAL_DIR_SIG;
        eocd.num_entries_disk = static_cast<uint16_t>(entries_.size());
        eocd.num_entries_total = static_cast<uint16_t>(entries_.size());
        eocd.central_dir_size = central_dir_size;
        eocd.central_dir_offset = central_dir_offset;
        
        file_.write(reinterpret_cast<const char*>(&eocd), sizeof(eocd));
        file_.close();
      }
      
    private:
      struct EntryRecord
      {
        std::string filename;
        uint32_t offset;
        uint32_t size;
        uint32_t crc32;
      };
      
      std::ofstream file_;
      std::vector<EntryRecord> entries_;
    };

  } // namespace Internal

  // ============================================================================
  // MzPeakFile implementation
  // ============================================================================

  MzPeakFile::MzPeakFile() = default;

  MzPeakFile::~MzPeakFile() = default;

  PeakFileOptions& MzPeakFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& MzPeakFile::getOptions() const
  {
    return options_;
  }

  void MzPeakFile::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  void MzPeakFile::setConfig(const MzPeakConfig& config)
  {
    config_ = config;
  }

  const MzPeakFile::MzPeakConfig& MzPeakFile::getConfig() const
  {
    return config_;
  }

  void MzPeakFile::load(const String& filename, MapType& map)
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    map.reset();
    
    // Set log if logging is enabled
    startProgress(0, 3, "Loading mzPeak file");
    
    // Step 1: Read file index
    setProgress(0);
    auto entries = Internal::readZipIndex(filename);
    
    // Find required files in archive
    bool has_spectrum_metadata = false;
    bool has_spectrum_data = false;
    
    for (const auto& entry : entries)
    {
      if (entry.filename == "spectra_metadata.mzpeak")
      {
        has_spectrum_metadata = true;
      }
      else if (entry.filename == "spectra_data.mzpeak")
      {
        has_spectrum_data = true;
      }
    }
    
    if (!has_spectrum_metadata || !has_spectrum_data)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
        "Required mzPeak files not found in archive (spectra_metadata.mzpeak, spectra_data.mzpeak)");
    }
    
    // Step 2: Read spectrum metadata
    setProgress(1);
    readSpectrumMetadata_(filename, map);
    
    // Step 3: Read spectrum data
    setProgress(2);
    readSpectrumData_(filename, map);
    
    // Step 4: Read chromatograms if present
    setProgress(3);
    readChromatograms_(filename, map);
    
    endProgress();
    
    OPENMS_LOG_INFO << "Loaded " << map.getNrSpectra() << " spectra and " 
                    << map.getNrChromatograms() << " chromatograms from " << filename << std::endl;
  }

  void MzPeakFile::store(const String& filename, const MapType& map) const
  {
    startProgress(0, 3, "Storing mzPeak file");
    
    // Create ZIP writer
    Internal::ZipWriter writer(filename);
    
    // Step 1: Write spectrum metadata
    setProgress(1);
    writeSpectrumMetadata_(filename, map);
    
    // Step 2: Write spectrum data
    setProgress(2);
    writeSpectrumData_(filename, map);
    
    // Step 3: Write chromatograms if enabled and present
    setProgress(3);
    if (config_.write_chromatograms && !map.getChromatograms().empty())
    {
      writeChromatograms_(filename, map);
    }
    
    endProgress();
    
    OPENMS_LOG_INFO << "Stored " << map.getNrSpectra() << " spectra and " 
                    << map.getNrChromatograms() << " chromatograms to " << filename << std::endl;
  }

  void MzPeakFile::transform(const String& filename_in, 
                             Interfaces::IMSDataConsumer* consumer, 
                             bool /* skip_full_count */, 
                             bool /* skip_first_pass */)
  {
    // Load the full experiment
    MapType exp;
    load(filename_in, exp);
    
    // Set expected size
    consumer->setExpectedSize(exp.getNrSpectra(), exp.getNrChromatograms());
    consumer->setExperimentalSettings(exp);
    
    // Feed spectra to consumer
    for (auto& spec : exp)
    {
      consumer->consumeSpectrum(spec);
    }
    
    // Feed chromatograms to consumer
    for (auto& chrom : exp.getChromatograms())
    {
      consumer->consumeChromatogram(chrom);
    }
  }

  void MzPeakFile::readSpectrumMetadata_(const String& archive_path, MapType& map)
  {
    auto entries = Internal::readZipIndex(archive_path);
    
    for (const auto& entry : entries)
    {
      if (entry.filename == "spectra_metadata.mzpeak")
      {
        auto data = Internal::extractZipEntry(archive_path, entry);
        
        // Create Arrow buffer from data
        auto buffer = std::make_shared<arrow::Buffer>(data.data(), data.size());
        auto buffer_reader = std::make_shared<arrow::io::BufferReader>(buffer);
        
        // Open Parquet file
        std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
        PARQUET_ASSIGN_OR_THROW(arrow_reader,
          parquet::arrow::OpenFile(buffer_reader, arrow::default_memory_pool())
        );
        
        // Read entire table
        std::shared_ptr<arrow::Table> table;
        PARQUET_THROW_NOT_OK(arrow_reader->ReadTable(&table));
        
        // Get spectrum column data
        // Expected columns: spectrum.index, spectrum.id, spectrum.time, spectrum.ms_level, etc.
        auto schema = table->schema();
        
        // Pre-allocate spectra
        map.resize(table->num_rows());
        
        // Find the 'spectrum' struct column
        int spectrum_col_idx = schema->GetFieldIndex("spectrum");
        if (spectrum_col_idx < 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path,
            "Could not find 'spectrum' column in spectra_metadata.mzpeak");
        }
        
        auto spectrum_col = table->column(spectrum_col_idx);
        
        // Process each row/spectrum
        for (int64_t i = 0; i < table->num_rows(); ++i)
        {
          MSSpectrum& spec = map[i];
          
          // Extract metadata from nested struct
          // For now, set minimal metadata - expand based on actual column structure
          auto chunk = std::static_pointer_cast<arrow::StructArray>(spectrum_col->chunk(0));
          
          // Get nested columns
          auto index_arr = std::static_pointer_cast<arrow::UInt64Array>(chunk->GetFieldByName("index"));
          auto time_arr = std::static_pointer_cast<arrow::DoubleArray>(chunk->GetFieldByName("time"));
          auto ms_level_arr = std::static_pointer_cast<arrow::UInt8Array>(chunk->GetFieldByName("ms_level"));
          
          if (time_arr && !time_arr->IsNull(i))
          {
            spec.setRT(time_arr->Value(i) * 60.0); // mzPeak uses minutes, OpenMS uses seconds
          }
          
          if (ms_level_arr && !ms_level_arr->IsNull(i))
          {
            spec.setMSLevel(ms_level_arr->Value(i));
          }
          
          // Set native ID from index
          if (index_arr && !index_arr->IsNull(i))
          {
            spec.setNativeID("spectrum=" + std::to_string(index_arr->Value(i)));
          }
        }
        
        break;
      }
    }
  }

  void MzPeakFile::readSpectrumData_(const String& archive_path, MapType& map)
  {
    auto entries = Internal::readZipIndex(archive_path);
    
    for (const auto& entry : entries)
    {
      if (entry.filename == "spectra_data.mzpeak")
      {
        auto data = Internal::extractZipEntry(archive_path, entry);
        
        // Create Arrow buffer from data
        auto buffer = std::make_shared<arrow::Buffer>(data.data(), data.size());
        auto buffer_reader = std::make_shared<arrow::io::BufferReader>(buffer);
        
        // Open Parquet file
        std::unique_ptr<parquet::arrow::FileReader> arrow_reader;
        PARQUET_ASSIGN_OR_THROW(arrow_reader,
          parquet::arrow::OpenFile(buffer_reader, arrow::default_memory_pool())
        );
        
        // Read entire table
        std::shared_ptr<arrow::Table> table;
        PARQUET_THROW_NOT_OK(arrow_reader->ReadTable(&table));
        
        auto schema = table->schema();
        
        // Find columns for spectrum index, m/z, and intensity
        int idx_col = schema->GetFieldIndex("spectrum_index");
        int mz_col = schema->GetFieldIndex("mz");
        int intensity_col = schema->GetFieldIndex("intensity");
        
        if (idx_col < 0 || mz_col < 0 || intensity_col < 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive_path,
            "Could not find required columns (spectrum_index, mz, intensity) in spectra_data.mzpeak");
        }
        
        auto idx_chunks = table->column(idx_col);
        auto mz_chunks = table->column(mz_col);
        auto intensity_chunks = table->column(intensity_col);
        
        // Process data in point layout format
        // Each row is a single peak: (spectrum_index, mz, intensity)
        int64_t current_idx = -1;
        std::vector<double> mz_values;
        std::vector<double> intensity_values;
        
        for (int chunk_idx = 0; chunk_idx < idx_chunks->num_chunks(); ++chunk_idx)
        {
          auto idx_arr = std::static_pointer_cast<arrow::UInt64Array>(idx_chunks->chunk(chunk_idx));
          
          // Handle both float32 and float64 for m/z
          std::shared_ptr<arrow::DoubleArray> mz_arr_f64;
          std::shared_ptr<arrow::FloatArray> mz_arr_f32;
          if (mz_chunks->chunk(chunk_idx)->type_id() == arrow::Type::DOUBLE)
          {
            mz_arr_f64 = std::static_pointer_cast<arrow::DoubleArray>(mz_chunks->chunk(chunk_idx));
          }
          else
          {
            mz_arr_f32 = std::static_pointer_cast<arrow::FloatArray>(mz_chunks->chunk(chunk_idx));
          }
          
          // Handle both float32 and float64 for intensity
          std::shared_ptr<arrow::DoubleArray> int_arr_f64;
          std::shared_ptr<arrow::FloatArray> int_arr_f32;
          if (intensity_chunks->chunk(chunk_idx)->type_id() == arrow::Type::DOUBLE)
          {
            int_arr_f64 = std::static_pointer_cast<arrow::DoubleArray>(intensity_chunks->chunk(chunk_idx));
          }
          else
          {
            int_arr_f32 = std::static_pointer_cast<arrow::FloatArray>(intensity_chunks->chunk(chunk_idx));
          }
          
          for (int64_t i = 0; i < idx_arr->length(); ++i)
          {
            int64_t spec_idx = idx_arr->Value(i);
            
            // If we moved to a new spectrum, save the previous one
            if (spec_idx != current_idx && current_idx >= 0)
            {
              if (static_cast<Size>(current_idx) < map.size())
              {
                auto& spec = map[current_idx];
                spec.resize(mz_values.size());
                for (Size j = 0; j < mz_values.size(); ++j)
                {
                  spec[j].setMZ(mz_values[j]);
                  spec[j].setIntensity(static_cast<Peak1D::IntensityType>(intensity_values[j]));
                }
              }
              mz_values.clear();
              intensity_values.clear();
            }
            
            current_idx = spec_idx;
            
            // Get m/z value
            double mz = mz_arr_f64 ? mz_arr_f64->Value(i) : static_cast<double>(mz_arr_f32->Value(i));
            double intensity = int_arr_f64 ? int_arr_f64->Value(i) : static_cast<double>(int_arr_f32->Value(i));
            
            mz_values.push_back(mz);
            intensity_values.push_back(intensity);
          }
        }
        
        // Don't forget the last spectrum
        if (current_idx >= 0 && static_cast<Size>(current_idx) < map.size())
        {
          auto& spec = map[current_idx];
          spec.resize(mz_values.size());
          for (Size j = 0; j < mz_values.size(); ++j)
          {
            spec[j].setMZ(mz_values[j]);
            spec[j].setIntensity(static_cast<Peak1D::IntensityType>(intensity_values[j]));
          }
        }
        
        break;
      }
    }
  }

  void MzPeakFile::readChromatograms_(const String& archive_path, MapType& map)
  {
    auto entries = Internal::readZipIndex(archive_path);
    
    bool has_chrom_metadata = false;
    bool has_chrom_data = false;
    
    for (const auto& entry : entries)
    {
      if (entry.filename == "chromatograms_metadata.mzpeak")
      {
        has_chrom_metadata = true;
      }
      else if (entry.filename == "chromatograms_data.mzpeak")
      {
        has_chrom_data = true;
      }
    }
    
    if (!has_chrom_metadata || !has_chrom_data)
    {
      // Chromatograms are optional
      return;
    }
    
    // TODO: Implement chromatogram reading similar to spectrum reading
    OPENMS_LOG_WARN << "Chromatogram reading not yet fully implemented for mzPeak format" << std::endl;
  }

  void MzPeakFile::writeSpectrumMetadata_(const String& /* archive_path */, const MapType& /* map */) const
  {
    // TODO: Implement spectrum metadata writing
    // This requires building Arrow/Parquet tables with the correct schema
    OPENMS_LOG_WARN << "Spectrum metadata writing not yet fully implemented for mzPeak format" << std::endl;
  }

  void MzPeakFile::writeSpectrumData_(const String& /* archive_path */, const MapType& /* map */) const
  {
    // TODO: Implement spectrum data writing
    // This requires building Arrow/Parquet tables with the correct schema
    OPENMS_LOG_WARN << "Spectrum data writing not yet fully implemented for mzPeak format" << std::endl;
  }

  void MzPeakFile::writeChromatograms_(const String& /* archive_path */, const MapType& /* map */) const
  {
    // TODO: Implement chromatogram writing
    OPENMS_LOG_WARN << "Chromatogram writing not yet fully implemented for mzPeak format" << std::endl;
  }

} // namespace OpenMS

#endif // WITH_PARQUET

