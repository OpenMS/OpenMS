// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cstdint>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS
{
  class MSSpectrum;

  /**
    @brief Dataset-level metadata for imzML mass spectrometry imaging files.

    imzML (HUPO-PSI standard v1.1.0) extends mzML with a companion binary
    file (.ibd) and IMS ontology CV terms (IMS:*) that describe the imaging
    geometry, binary layout, and file checksum.

    @p ImzMLFile::load stores these fields on the loaded @p MSExperiment as
    MetaValues (see @p ImzMLFile documentation).

    @see ImzMLFile, OnDiscImzMLExperiment, Internal::ImzMLHandler
    @ingroup FileIO
  */
  struct OPENMS_DLLAPI ImzMLMeta
  {
    // -------------------------------------------------------------------------
    // Imaging geometry
    // -------------------------------------------------------------------------

    uint32_t max_count_x  {0};  ///< Max pixel column (IMS:1000042)
    uint32_t max_count_y  {0};  ///< Max pixel row    (IMS:1000043)
    uint32_t max_count_z  {1};  ///< Max depth slice  (usually 1 for 2-D datasets)
    double   pixel_size_x {-1}; ///< Physical pixel width  in @p pixel_size_unit (IMS:1000046), -1 if unknown
    double   pixel_size_y {-1}; ///< Physical pixel height in @p pixel_size_unit (IMS:1000047), -1 if unknown
    std::string pixel_size_unit {"unknown"}; ///< UO length unit name of the pixel size, "unknown" if no pixel size
    double   max_dim_x    {0};  ///< Physical x extent in µm (IMS:1000044)
    double   max_dim_y    {0};  ///< Physical y extent in µm (IMS:1000045)

    // -------------------------------------------------------------------------
    // Imaging mode
    // -------------------------------------------------------------------------

    /// "continuous" (shared m/z array) or "processed" (per-spectrum m/z).
    std::string imaging_mode;

    // -------------------------------------------------------------------------
    // IBD companion file
    // -------------------------------------------------------------------------

    std::string ibd_file_path; ///< Absolute path to the .ibd file
    std::string ibd_sha1;      ///< SHA-1 checksum (IMS:1000091), empty if absent
    std::string ibd_md5;       ///< MD5  checksum (IMS:1000090), empty if absent
    std::string uuid;          ///< Dataset UUID (IMS:1000080)

    // -------------------------------------------------------------------------
    // Array data types (first occurrence, dataset-level summary)
    // -------------------------------------------------------------------------

    std::string mz_data_type;   ///< "float32" | "float64" | "int32" | "int64"
    std::string int_data_type;  ///< same set as mz_data_type

    // -------------------------------------------------------------------------
    // Acquisition geometry
    // -------------------------------------------------------------------------

    std::string scan_pattern;        ///< "top down" | "bottom up" (IMS:1000401/402)
    std::string scan_direction;      ///< "flyback" | "meander" | "horizontal" | "vertical"
    std::string line_scan_direction; ///< "left-right" | "right-left" (IMS:1000491/492)
    std::string polarity;            ///< "positive" | "negative"  (MS:1000130/129)
  };

  /**
    @brief SI exponent of the UO length unit @p unit_accession, i.e. its power of ten relative to meter.

    The imzML pixel size (IMS:1000046/47) may use a UO length unit, while the max
    dimension (IMS:1000044/45) is fixed to micrometer. Only meter, centimeter, millimeter,
    micrometer, nanometer, angstrom and picometer are known; this also serves as the
    check which pixel size units are supported.

    @param[in] unit_accession UO accession of the length unit, e.g. "UO:0000018" (nanometer)
    @return SI exponent, e.g. -9 for nanometer, or std::nullopt if the unit is not known
  */
  OPENMS_DLLAPI std::optional<int> lengthUnitExponent(const std::string& unit_accession);


  /**
    @brief Per-spectrum binary index entry for an imzML dataset.

    Records the byte offset and element count for the m/z and intensity
    arrays in the companion .ibd file, together with the pixel (x,y,z)
    coordinates and zlib-compression flags (@c MS:1000574).  Built during
    ImzMLHandler parsing and stored inside OnDiscImzMLExperiment so that
    random-access spectrum decoding requires only a single fseek + fread
    per call. Compressed m/z, intensity, and auxiliary arrays (any child of
    @c MS:1000572 other than uncompressed @c MS:1000576) are rejected at
    decode time (inflate / numpress of external arrays is not supported).

    Optional @p aux entries describe additional external arrays (ion mobility
    and other @c FloatDataArray values). Both in-memory @p ImzMLFile::load and
    @p OnDiscImzMLExperiment::getSpectrum() attach those as @c FloatDataArray
    on the returned spectrum when present.

    @ingroup FileIO
  */
  struct OPENMS_DLLAPI ImzMLSpectrumIndex
  {
    /// Scalar type identifier for binary array elements.
    enum class DataType : uint8_t { FLOAT32, FLOAT64, INT32, INT64, UNKNOWN };

    /**
      @brief Index entry for one auxiliary (non-m/z, non-intensity) external array.

      Used for ion mobility and other per-peak FloatDataArrays so on-disc and
      in-memory loaders expose the same viewer contract (@c containsIMData(),
      equal length to peaks).
    */
    struct AuxArray
    {
      std::string name;                 ///< Ontology / free-text array name (FloatDataArray::getName)
      std::string accession;            ///< MS accession (e.g. MS:1003006), may be empty
      std::string unit_accession;       ///< unitAccession on the array-identity cvParam (may be empty)
      uint64_t offset {0};              ///< IMS:1000102
      uint64_t length {0};              ///< IMS:1000103 element count
      uint64_t encoded_bytes {0};       ///< IMS:1000104 byte length on disk
      DataType type {DataType::UNKNOWN};
      bool compressed {false};          ///< Child of MS:1000572 other than MS:1000576 (unsupported)
    };

    int32_t  index      {0};               ///< 0-based document order
    uint32_t x          {0};               ///< Pixel column (1-based, IMS:1000050)
    uint32_t y          {0};               ///< Pixel row    (1-based, IMS:1000051)
    uint32_t z          {1};               ///< Depth slice  (1-based, IMS:1000052)
    uint64_t mz_offset  {0};              ///< Byte offset of m/z array (IMS:1000102)
    uint64_t mz_length  {0};              ///< Element count             (IMS:1000103)
    DataType mz_type    {DataType::UNKNOWN};
    bool mz_compressed {false};           ///< Compressed (not MS:1000576) on m/z array
    uint64_t int_offset {0};              ///< Byte offset of intensity array
    uint64_t int_length {0};              ///< Element count
    DataType int_type   {DataType::UNKNOWN};
    bool int_compressed {false};          ///< Compressed (not MS:1000576) on intensity array
    std::vector<AuxArray> aux;            ///< Extra external arrays (IM, …) in document order
  };

  /**
    @brief Shared .ibd binary array decode utilities for imzML handlers.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ImzMLBinaryIO
  {
  public:
    /**
      @brief Read an m/z array from the companion .ibd file.

      @param[in] ibd       Open binary file handle.
      @param[in] offset    Byte offset (IMS:1000102).
      @param[in] count     Element count (IMS:1000103).
      @param[in] dt        Scalar type of stored values.
      @param[out] out       Decoded m/z values.
      @param[in] ibd_path  Path used in error messages.

      @throws Exception::ParseError if seek or read fails.
    */
    static void readMzArray(FILE* ibd,
                            uint64_t offset,
                            uint64_t count,
                            ImzMLSpectrumIndex::DataType dt,
                            std::vector<double>& out,
                            const std::string& ibd_path);

    /**
      @brief Read an intensity array from the companion .ibd file.

      @throws Exception::ParseError if seek or read fails.
    */
    static void readIntArray(FILE* ibd,
                             uint64_t offset,
                             uint64_t count,
                             ImzMLSpectrumIndex::DataType dt,
                             std::vector<float>& out,
                             const std::string& ibd_path);

    /**
      @brief Read an auxiliary (non-m/z, non-intensity) array from the .ibd file.

      imzML permits additional @c binaryDataArray entries per spectrum beyond the
      mandatory m/z and intensity pair — most notably ion mobility arrays such as
      @c MS:1003006 "mean inverse reduced ion mobility array". Values are returned
      as @c float because OpenMS stores auxiliary spectrum data in
      @p MSSpectrum::FloatDataArray.

      @param[in] ibd        Open binary file handle.
      @param[in] offset     Byte offset (IMS:1000102).
      @param[in] count      Element count (IMS:1000103).
      @param[in] dt         Scalar type of stored values.
      @param[out] out       Decoded values as float (FloatDataArray storage).
      @param[in] ibd_path   Path used in error messages.
      @param[in] array_name Array name used in error messages (e.g. the CV term name).

      @note Compressed external arrays (any child of MS:1000572 other than
      uncompressed MS:1000576) are not decoded; callers must reject them before
      invoking this method. @c IMS:1000104 is then required to know how many
      compressed bytes to read.

      @throws Exception::ParseError if seek or read fails.
    */
    static void readAuxArray(FILE* ibd,
                             uint64_t offset,
                             uint64_t count,
                             ImzMLSpectrumIndex::DataType dt,
                             std::vector<float>& out,
                             const std::string& ibd_path,
                             const std::string& array_name);

    /**
      @brief Write a float32 array to the companion .ibd file (little-endian).

      @throws Exception::ParseError if write fails.
    */
    static void writeFloat32Array(FILE* ibd,
                                  const float* data,
                                  uint64_t count,
                                  const std::string& ibd_path);

    /**
      @brief Write m/z values as float32 to the companion .ibd file.

      @throws Exception::ParseError if write fails.
    */
    static void writeMzAsFloat32(FILE* ibd,
                                 const std::vector<double>& mz,
                                 const std::string& ibd_path);

    /**
      @brief Write a float64 array to the companion .ibd file (little-endian).

      @throws Exception::ParseError if write fails.
    */
    static void writeFloat64Array(FILE* ibd,
                                  const double* data,
                                  uint64_t count,
                                  const std::string& ibd_path);

    /**
      @brief Write m/z values as float64 to the companion .ibd file.

      @throws Exception::ParseError if write fails.
    */
    static void writeMzAsFloat64(FILE* ibd,
                                 const std::vector<double>& mz,
                                 const std::string& ibd_path);
  };

  /**
    @brief Per-spectrum rules shared by every imzML access path (in-memory, streaming,
    metadata-only, index-only and @p OnDiscImzMLExperiment), so they cannot drift apart.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ImzMLSpectrumRules
  {
  public:
    /**
      @brief Acquisition (scan start) time of an imaging spectrum in seconds.

      Imaging spectra keep RT at -1 and store the time as the MetaValue
      @p imzml:scan_start_time. An RT the caller set (>= 0) wins over the MetaValue;
      a non-numeric or negative MetaValue counts as no time (a scan start time is
      relative to the start of the run, so a negative one is a placeholder).

      @return The time (>= 0), or the (negative) RT if the spectrum has none.
    */
    static double scanStartTime(const MSSpectrum& s);

    /**
      @brief Whether an external auxiliary array is read from the .ibd for a spectrum
      with @p n_peaks peaks (it needs an identity, a non-zero length equal to the peak
      count and a known binary data type). The one rule for every access path.

      With @p warn, a skipped array (other than one of length 0) is logged with the pixel
      @p x / @p y / @p z.
    */
    static bool isLoadableAuxArray(const std::string& name, uint64_t length,
                                   ImzMLSpectrumIndex::DataType type, uint64_t n_peaks,
                                   bool warn = false, uint32_t x = 0, uint32_t y = 0, uint32_t z = 0);

    /**
      @brief Set the IM peak type after the auxiliary arrays are known: with an ion
      mobility array, @c UNKNOWN becomes @c IM_PROFILE (a type the file declares, e.g.
      @c MS:1003441, is kept); without one, the type is @c UNKNOWN.
    */
    static void applyIMPeakType(MSSpectrum& s, bool has_im_array);

    /**
      @brief Attach the values of an external auxiliary array read from the .ibd.

      They go into the first still-empty FloatDataArray named @p name (MzMLHandler creates
      one per auxiliary binaryDataArray, with its metadata), so repeated names map 1:1 in
      order; without one, a new array named @p name (and @p unit_accession) is added.
    */
    static void attachAuxArray(MSSpectrum& s, const std::string& name, const std::string& unit_accession,
                               const std::vector<float>& values);

    /// Drop the empty Float/Integer/StringDataArrays left by MzMLHandler for aux arrays that were not filled.
    static void pruneEmptyDataArrays(MSSpectrum& s);
  };

} // namespace OpenMS
