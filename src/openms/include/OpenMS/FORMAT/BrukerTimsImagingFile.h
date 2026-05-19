// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

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
  */
  class OPENMS_DLLAPI BrukerTimsImagingFile : public ProgressLogger
  {
  public:
    /// One row of MaldiFrameInfo (2D — z is intentionally not exposed).
    struct MaldiPixel
    {
      Int frame_id = 0;
      Int x = 0;
      Int y = 0;
    };

    /// Processing and export configuration.
    struct Config
    {
      /// Configuration for the inner BrukerTimsFile (calibration, centroiding,
      /// etc.). @c export_mode is forced to FRAME and @c load_ms1 to true
      /// before delegation.
      BrukerTimsFile::Config inner_config;

      /// Reject datasets whose MaldiApplicationType is not "Imaging".
      bool strict_imaging_only = true;
    };

    BrukerTimsImagingFile() = default;

    /// Load a MALDI imaging .d folder into @p exp.
    /// @throws Exception::FileNotReadable if the .d folder or analysis.tdf is missing.
    /// @throws Exception::InvalidValue if the dataset is not MALDI imaging and strict_imaging_only is true,
    ///         or if the dataset contains multiple Z slices (out of scope; load each separately).
    /// @throws Exception::ParseError if MaldiFrameInfo is absent or empty.
    void load(const String& path, MSImagingExperiment& exp);
    /// @overload
    void load(const String& path, MSImagingExperiment& exp, const Config& config);

    /// Quick probe: returns true iff GlobalMetadata.MaldiApplicationType == "Imaging".
    /// Returns false on any error (missing file, unreadable, etc.) so callers can
    /// use it as a non-throwing duck-type check.
    static bool isImagingDataset(const String& path);

    /// Reads @c MaldiFrameInfo and returns rows in @c Frame.Id order.
    /// @throws Exception::FileNotReadable if @c analysis.tdf cannot be opened.
    /// @throws Exception::ParseError if @c MaldiFrameInfo is missing or empty.
    static std::vector<MaldiPixel> readMaldiFrameInfo(const String& d_folder);

    /// Reads a single @c GlobalMetadata key.
    /// @returns the value if present, otherwise an empty string.
    /// @throws Exception::FileNotReadable if @c analysis.tdf cannot be opened.
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
