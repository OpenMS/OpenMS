// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <string>
#include <vector>

namespace OpenMS
{

  class MSImagingExperiment;
  class MSImagingGeometry;

  /**
    @brief File adapter for imzML 1.1.0 mass spectrometry imaging data.

    imzML extends mzML with a companion binary file and IMS ontology terms.
    Data are split across two files on disk:

    - @c .imzML — XML metadata (mzML 1.1.0 schema + IMS CV terms)
    - @c .ibd   — external m/z and intensity arrays

    The @c .ibd path is derived from the @c .imzML path (same basename,
    extension replaced with @c .ibd, case-insensitive).

    Parsing reuses @p MzMLHandler via @p Internal::ImzMLHandler, which
    decodes IMS-specific terms and reads peak data from the @c .ibd file.

    @b Loading modes
    - In-memory: @p load(filename, MSImagingExperiment&) decodes the spectra and
      builds @p MSImagingGeometry for random access by (x, y) in RAM (zero-based
      coordinates; imzML file coordinates are 1-based). Use
      @p MSImagingExperiment::getMSExperiment() for the flat spectra.
    - Streaming: @p load(filename, IMSDataConsumer&) forwards decoded spectra to
      the consumer without retaining them in @p ImzMLFile. Delivery is batched
      after the @c spectrumList XML section is parsed (same as @p MzMLHandler),
      not on a per-startElement basis.
    - Index only: @p loadSpectraIndex() builds a byte-offset index for
      @p OnDiscImzMLExperiment.

    @b Imaging metadata on @p MSExperiment
    After @p load(), dataset fields are stored as MetaValues, e.g.
    @p imzml:imaging_mode, @p imzml:ibd_path, @p imzml:uuid, checksums,
    @p imzml:max_count_x/y/z, pixel sizes, scan geometry, polarity, and array
    data types when present in the file.
    Per-spectrum pixel coordinates use @p imzml:x, @p imzml:y, @p imzml:z.

    @ingroup FileIO

    @see MzMLFile, OnDiscImzMLExperiment, MSImagingExperiment, ImzMLMeta, ImzMLSpectrumIndex
  */
  class OPENMS_DLLAPI ImzMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
  public:

    /** @name Constructors and destructor */
    //@{
    ImzMLFile();
    ~ImzMLFile() override = default;
    //@}

    /** @name Options */
    //@{
    PeakFileOptions& getOptions();
    const PeakFileOptions& getOptions() const;
    void setOptions(const PeakFileOptions& options);
    //@}

    /** @name Load */
    //@{
    /**
      @brief Load an imzML file into an MSImagingExperiment for in-memory random access.

      Performs a full load, then builds @p MSImagingGeometry from @p imzml:x/y
      MetaValues on each spectrum (imzML 1-based coordinates are mapped to
      zero-based pixel indices). Only spectra with @p imzml:z = 1 are indexed
      (2-D datasets).

      Use @p MSImagingExperiment::getSpectrum(x, y) to fetch spectra by pixel
      without re-reading the @c .ibd file.

      @param[in] filename Path to the @c .imzML file.
      @param[out] exp       Target imaging experiment (replaced on success).

      @throws Exception::FileNotFound if the @c .imzML or @c .ibd cannot be opened.
      @throws Exception::ParseError   if the XML is malformed.
      @throws Exception::InvalidValue if pixel coordinates are invalid or duplicate.
    */
    void load(const std::string& filename, MSImagingExperiment& exp);

    /**
      @brief Build @p MSImagingGeometry from a loaded imzML @p MSExperiment.

      Reads dataset dimensions from @p imzml:max_count_x/y MetaValues when
      present and registers each spectrum that carries @p imzml:x/y (and
      @p imzml:z = 1). Coordinates are converted from imzML's 1-based convention
      to zero-based geometry indices.

      @param[in] exp  Experiment previously loaded from imzML (e.g. via @p load or @p FileHandler).
      @param[out] geom Geometry to populate (cleared first).

      @throws Exception::InvalidValue on invalid or duplicate pixel coordinates.
    */
    static void buildImagingGeometry(const MSExperiment& exp, MSImagingGeometry& geom);

    /**
      @brief Build @p MSImagingGeometry directly from a parsed imzML spectrum index.

      This is the source-of-truth path used by the loaders: it derives the pixel grid
      straight from the per-spectrum coordinates captured during parsing (no reliance on
      @p imzml:x/y MetaValues), mirroring how BrukerTimsImagingFile populates its geometry.
      Each index entry's 1-based (x, y) with z = 1 is registered as a zero-based pixel bound
      to that spectrum; dimensions and pixel size come from @p meta. Out-of-grid or < 1
      coordinates are warned and skipped; duplicates raise Exception::InvalidValue.

      @param[in]  index Per-spectrum index entries (e.g. from @p loadSpectraIndex or the loader).
      @param[in]  meta  Dataset metadata (grid dimensions, pixel size).
      @param[out] geom  Geometry to populate (cleared first).

      @throws Exception::InvalidValue on duplicate pixel coordinates.
    */
    static void buildImagingGeometry(const std::vector<ImzMLSpectrumIndex>& index,
                                     const ImzMLMeta& meta,
                                     MSImagingGeometry& geom);

    /**
      @brief Stream-load an imzML file into a consumer.

      Decoded spectra are passed to @p consumer.consumeSpectrum() after the
      @c spectrumList section has been parsed (batched delivery via
      @p MzMLHandler). Peak data are not retained in @p ImzMLFile itself.

      @param[in] filename Path to the @c .imzML file.
      @param[in] consumer Receiver for spectra and experimental settings.

      @throws Exception::FileNotFound if the @c .imzML or @c .ibd cannot be opened.
      @throws Exception::ParseError   if the XML is malformed.
    */
    void load(const std::string& filename, Interfaces::IMSDataConsumer& consumer);

    /**
      @brief Parse XML and build a per-spectrum @c .ibd index without loading peaks.

      Populates @p meta and @p index for use with @p OnDiscImzMLExperiment.
      Opens the companion @c .ibd and records its path in @p meta.ibd_file_path.

      @param[in] filename Path to the @c .imzML file.
      @param[out] meta     Dataset-level imaging metadata.
      @param[out] index    Per-spectrum index entries (offsets into @c .ibd).
      @param[in] ibd_path  Optional explicit path to the companion @c .ibd; inferred from
                           @p filename when empty.

      @throws Exception::FileNotFound if the @c .imzML or @c .ibd cannot be opened.
      @throws Exception::ParseError   if the XML is malformed.
    */
    void loadSpectraIndex(const std::string& filename,
                          ImzMLMeta& meta,
                          std::vector<ImzMLSpectrumIndex>& index,
                          const std::string& ibd_path = "");
    //@}

    /** @name Validation and storage */
    //@{
    /**
      @brief Validate an imzML file against the mzML 1.1.0 XML schema.

      @param[in] filename Path to the @c .imzML file.
      @param[out] os       Stream receiving validation messages.
      @return @c true if the file passes schema validation.
    */
    bool isValid(const std::string& filename, std::ostream& os);

    /**
      @brief Store an experiment as imzML (XML + companion .ibd).

      Writes external binary arrays (float32 or float64 via @p PeakFileOptions) with a
      16-byte UUID header in the @c .ibd file linked to IMS:1000080 in the XML. Continuous mode is selected when
      @p imzml:imaging_mode is @c continuous or all spectra share an identical
      m/z axis; otherwise processed mode is used.

      Each spectrum must carry @p imzml:x and @p imzml:y MetaValues (1-based imzML
      pixel coordinates). Dataset imaging metadata (@p imzml:scan_pattern,
      @p imzml:polarity, grid dimensions, etc.) and per-spectrum RT/MS level are
      written when present. @p PeakFileOptions control binary precision and
      optional spectrum/peak filtering during export.

      @param[in] filename Path to the output @c .imzML file.
      @param[in] exp      Experiment with spectra and optional imzML MetaValues.

      @throws Exception::MissingInformation if @p exp has no spectra or lacks imzml:x/y.
      @throws Exception::InvalidValue if pixel coordinates are invalid or duplicate.
      @throws Exception::InvalidParameter if continuous export is requested but spectra are incompatible.
      @throws Exception::UnableToCreateFile if output files cannot be written.
      @throws Exception::ParseError if binary array serialization fails.
    */
    void store(const std::string& filename, const MSExperiment& exp) const;

    /**
      @brief Store an MSImagingExperiment as imzML, taking pixel coordinates and grid
             dimensions from its MSImagingGeometry (the source of truth).

      Unlike the @p MSExperiment overload, this does not require @p imzml:x/y MetaValues on
      the spectra: the spatial information is read from @p exp.getGeometry(), so any
      MSImagingExperiment can be written — including ones built without those MetaValues
      (e.g. from @p BrukerTimsImagingFile). Dataset-level imzML metadata already present on
      the wrapped experiment (imaging mode, data types, scan geometry, ...) is preserved;
      grid dimensions and pixel size are taken from the geometry.

      @param[in] filename Path to the output @c .imzML file.
      @param[in] exp      Imaging experiment to store.

      @throws Exception::InvalidValue if a geometry pixel references a spectrum index outside
              the experiment, or pixel coordinates are invalid or duplicate.
      @throws Exception::UnableToCreateFile if output files cannot be written.
    */
    void store(const std::string& filename, const MSImagingExperiment& exp) const;
    //@}

  private:

    void loadImpl_(const std::string& filename,
                   Interfaces::IMSDataConsumer* consumer,
                   MSExperiment& meta_exp,
                   ImzMLMeta* out_meta,
                   std::vector<ImzMLSpectrumIndex>* out_index,
                   bool index_only,
                   const std::string& ibd_path_override = "");

    static std::string inferIbdPath_(const std::string& imzml_path);

    PeakFileOptions options_;
  };

} // namespace OpenMS
