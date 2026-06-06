// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdint>
#include <cstdio>
#include <vector>

namespace OpenMS
{

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
    double   pixel_size_x {0};  ///< Physical pixel width  in µm (IMS:1000046)
    double   pixel_size_y {0};  ///< Physical pixel height in µm (IMS:1000047)
    double   max_dim_x    {0};  ///< Physical x extent in µm (IMS:1000044)
    double   max_dim_y    {0};  ///< Physical y extent in µm (IMS:1000045)

    // -------------------------------------------------------------------------
    // Imaging mode
    // -------------------------------------------------------------------------

    /// "continuous" (shared m/z array) or "processed" (per-spectrum m/z).
    String imaging_mode;

    // -------------------------------------------------------------------------
    // IBD companion file
    // -------------------------------------------------------------------------

    String ibd_file_path; ///< Absolute path to the .ibd file
    String ibd_sha1;      ///< SHA-1 checksum (IMS:1000091), empty if absent
    String ibd_md5;       ///< MD5  checksum (IMS:1000090), empty if absent
    String uuid;          ///< Dataset UUID (IMS:1000080)

    // -------------------------------------------------------------------------
    // Array data types (first occurrence, dataset-level summary)
    // -------------------------------------------------------------------------

    String mz_data_type;   ///< "float32" | "float64" | "int32" | "int64"
    String int_data_type;  ///< same set as mz_data_type

    // -------------------------------------------------------------------------
    // Acquisition geometry
    // -------------------------------------------------------------------------

    String scan_pattern;        ///< "top down" | "bottom up" (IMS:1000401/402)
    String scan_direction;      ///< "flyback" | "meander" | "horizontal" | "vertical"
    String line_scan_direction; ///< "left-right" | "right-left" (IMS:1000491/492)
    String polarity;            ///< "positive" | "negative"  (MS:1000130/129)
  };


  /**
    @brief Per-spectrum binary index entry for an imzML dataset.

    Records the byte offset and element count for both the m/z and intensity
    arrays in the companion .ibd file, together with the pixel (x,y,z)
    coordinates.  Built during ImzMLHandler parsing and stored inside
    OnDiscImzMLExperiment so that random-access spectrum decoding requires
    only a single fseek + fread per call.

    @ingroup FileIO
  */
  struct OPENMS_DLLAPI ImzMLSpectrumIndex
  {
    /// Scalar type identifier for binary array elements.
    enum class DataType : uint8_t { FLOAT32, FLOAT64, INT32, INT64, UNKNOWN };

    int32_t  index      {0};               ///< 0-based document order
    uint32_t x          {0};               ///< Pixel column (1-based, IMS:1000050)
    uint32_t y          {0};               ///< Pixel row    (1-based, IMS:1000051)
    uint32_t z          {1};               ///< Depth slice  (1-based, IMS:1000052)
    uint64_t mz_offset  {0};              ///< Byte offset of m/z array (IMS:1000102)
    uint64_t mz_length  {0};              ///< Element count             (IMS:1000103)
    DataType mz_type    {DataType::UNKNOWN};
    uint64_t int_offset {0};              ///< Byte offset of intensity array
    uint64_t int_length {0};              ///< Element count
    DataType int_type   {DataType::UNKNOWN};
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
                            const String& ibd_path);

    /**
      @brief Read an intensity array from the companion .ibd file.

      @throws Exception::ParseError if seek or read fails.
    */
    static void readIntArray(FILE* ibd,
                             uint64_t offset,
                             uint64_t count,
                             ImzMLSpectrumIndex::DataType dt,
                             std::vector<float>& out,
                             const String& ibd_path);

    /**
      @brief Write a float32 array to the companion .ibd file (little-endian).

      @throws Exception::ParseError if write fails.
    */
    static void writeFloat32Array(FILE* ibd,
                                  const float* data,
                                  uint64_t count,
                                  const String& ibd_path);

    /**
      @brief Write m/z values as float32 to the companion .ibd file.

      @throws Exception::ParseError if write fails.
    */
    static void writeMzAsFloat32(FILE* ibd,
                                 const std::vector<double>& mz,
                                 const String& ibd_path);

    /**
      @brief Write a float64 array to the companion .ibd file (little-endian).

      @throws Exception::ParseError if write fails.
    */
    static void writeFloat64Array(FILE* ibd,
                                  const double* data,
                                  uint64_t count,
                                  const String& ibd_path);

    /**
      @brief Write m/z values as float64 to the companion .ibd file.

      @throws Exception::ParseError if write fails.
    */
    static void writeMzAsFloat64(FILE* ibd,
                                 const std::vector<double>& mz,
                                 const String& ibd_path);
  };

} // namespace OpenMS
