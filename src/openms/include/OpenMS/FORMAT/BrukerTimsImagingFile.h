// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Reader for Bruker TimsTOF MALDI imaging .d directories.

    Reads pixel coordinates from the @c MaldiFrameInfo SQL table inside
    @c analysis.tdf and binds them to per-frame spectra produced by
    BrukerTimsFile in FRAME export mode. No proprietary MALDI SDK calls
    are required — XY positions live entirely in the SQLite metadata.

    Coordinate convention: imzML / Bruker @c XIndexPos and @c YIndexPos
    are normalized to zero-based pixel coordinates inside MSImagingGeometry
    by subtracting the minimum observed index on each axis.

    MSImagingGeometry is strictly 2D. Multi-section MALDI datasets
    (multiple distinct @c ZIndexPos values) are rejected; load each
    section into its own MSImagingExperiment instead.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI BrukerTimsImagingFile : public ProgressLogger
  {
  public:
    /// One row of MaldiFrameInfo (2D — z is intentionally not exposed).
    struct MaldiPixel
    {
      Int frame_id = 0; ///< Bruker frame ID (1-based in .tdf).
      Int x = 0;        ///< Raw XIndexPos as stored in MaldiFrameInfo.
      Int y = 0;        ///< Raw YIndexPos as stored in MaldiFrameInfo.
    };

    /// Processing and export configuration.
    struct Config
    {
      /// Configuration for the inner BrukerTimsFile (calibration, centroiding,
      /// etc.). @c export_mode is forced to FRAME and @c load_ms1 to true
      /// before delegation.
      BrukerTimsFile::Config inner_config;

      /// Reject datasets whose @c MaldiApplicationType is not "Imaging".
      bool strict_imaging_only = true;
    };

    /// @brief Default constructor.
    BrukerTimsImagingFile() = default;

    /**
      @brief Loads a MALDI imaging .d folder into @p exp using default Config.
      @param[in]  path The path to the Bruker .d directory.
      @param[out] exp  Destination experiment; replaced (not appended).
      @throws Exception::FileNotReadable if the .d folder or @c analysis.tdf is missing.
      @throws Exception::InvalidValue if the dataset is not MALDI imaging and
              @c strict_imaging_only is true, or if the dataset contains
              multiple distinct @c ZIndexPos values (out of scope for the 2D
              MSImagingGeometry; load each section separately).
      @throws Exception::ParseError if @c MaldiFrameInfo is absent or empty.
    */
    void load(const String& path, MSImagingExperiment& exp);

    /**
      @brief Loads a MALDI imaging .d folder into @p exp.
      @param[in]  path   The path to the Bruker .d directory.
      @param[out] exp    Destination experiment; replaced (not appended).
      @param[in]  config Loader configuration; @c export_mode is forced to FRAME
                         and @c load_ms1 to true internally.
      @throws Exception::FileNotReadable if the .d folder or @c analysis.tdf is missing.
      @throws Exception::InvalidValue if the dataset is not MALDI imaging and
              @c strict_imaging_only is true, or if the dataset contains
              multiple distinct @c ZIndexPos values.
      @throws Exception::ParseError if @c MaldiFrameInfo is absent or empty.
    */
    void load(const String& path, MSImagingExperiment& exp, const Config& config);

    /**
      @brief Cheap, non-throwing probe for MALDI imaging mode.

      Returns false on any I/O error (missing .d folder, unreadable
      @c analysis.tdf, missing @c GlobalMetadata table) so callers can use this
      as a duck-type check. To distinguish "not imaging" from "unreadable",
      call @c readGlobalMetadataValue() instead, which throws on I/O failure.

      @param[in] path The path to the Bruker .d directory.
      @return true iff @c GlobalMetadata.MaldiApplicationType == "Imaging".
    */
    static bool isImagingDataset(const String& path);

    /**
      @brief Reads @c MaldiFrameInfo and returns rows in @c Frame.Id order.

      The @c ZIndexPos column, if present, is intentionally not read — the
      2D geometry contract is enforced by @c load(), which rejects multi-Z
      datasets up front.

      @param[in] d_folder The path to the Bruker .d directory.
      @return Pixel rows sorted ascending by frame id.
      @throws Exception::FileNotReadable if @c analysis.tdf cannot be opened.
      @throws Exception::ParseError if @c MaldiFrameInfo is missing or empty.
    */
    static std::vector<MaldiPixel> readMaldiFrameInfo(const String& d_folder);

    /**
      @brief Reads a single @c GlobalMetadata key.
      @param[in] d_folder The path to the Bruker .d directory.
      @param[in] key      Key to look up in the @c GlobalMetadata table.
      @return The stored value, or an empty string if the key is absent or
              the @c GlobalMetadata table does not exist.
      @throws Exception::FileNotReadable if @c analysis.tdf cannot be opened.
    */
    static String readGlobalMetadataValue(const String& d_folder, const String& key);

  private:
    /// Computes the 2D bounding box of the observed pixel coordinates.
    static void computeBoundingBox_(const std::vector<MaldiPixel>& rows,
                                    Int& x_min, Int& y_min,
                                    UInt& width, UInt& height);

    /// Resolves @c analysis.tdf inside a .d folder. Throws FileNotReadable if absent.
    static String resolveTdfPath_(const String& d_folder);
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
