// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <memory>

namespace OpenMS
{
  class MSChromatogramParquetConsumerImpl;

  /**
    @brief Writes chromatograms to a Parquet file with a PyProphet-compatible schema.

    The schema includes precursor/transition metadata, RT/intensity arrays and
    compression flags. Additional columns are run_id, source_file, and ms_level.

    The Parquet output has the following columns (one row per chromatogram):

      <table>
        <tr> <th BGCOLOR="#EBEBEB">Column</th> <th BGCOLOR="#EBEBEB">Type</th> <th BGCOLOR="#EBEBEB">Description</th> </tr>
        <tr> <td>RUN_ID</td> <td>int64</td> <td>Run identifier</td> </tr>
        <tr> <td>SOURCE_FILE</td> <td>string</td> <td>Input source filename</td> </tr>
        <tr> <td>MS_LEVEL</td> <td>int64</td> <td>MS level (1 for precursor traces, 2 for fragment traces)</td> </tr>
        <tr> <td>PRECURSOR_ID</td> <td>int64 (nullable)</td> <td>Precursor id</td> </tr>
        <tr> <td>TRANSITION_ID</td> <td>int64 (nullable)</td> <td>Transition id</td> </tr>
        <tr> <td>MODIFIED_SEQUENCE</td> <td>string (nullable)</td> <td>Modified peptide sequence</td> </tr>
        <tr> <td>PRECURSOR_CHARGE</td> <td>int64 (nullable)</td> <td>Precursor charge</td> </tr>
        <tr> <td>PRODUCT_CHARGE</td> <td>int64 (nullable)</td> <td>Product charge</td> </tr>
        <tr> <td>DETECTING_TRANSITION</td> <td>int64 (nullable)</td> <td>Detecting transition flag</td> </tr>
        <tr> <td>PRECURSOR_DECOY</td> <td>int64 (nullable)</td> <td>Precursor decoy flag</td> </tr>
        <tr> <td>PRODUCT_DECOY</td> <td>int64 (nullable)</td> <td>Product decoy flag</td> </tr>
        <tr> <td>TRANSITION_ORDINAL</td> <td>int64 (nullable)</td> <td>Transition ordinal</td> </tr>
        <tr> <td>TRANSITION_TYPE</td> <td>string (nullable)</td> <td>Transition type (e.g., y, b)</td> </tr>
        <tr> <td>ANNOTATION</td> <td>string (nullable)</td> <td>Transition annotation (e.g., y3^1)</td> </tr>
        <tr> <td>RT_DATA</td> <td>binary</td> <td>Compressed RT array</td> </tr>
        <tr> <td>INTENSITY_DATA</td> <td>binary</td> <td>Compressed intensity array</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>int64</td> <td>RT compression scheme id</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>int64</td> <td>Intensity compression scheme id</td> </tr>
      </table>

    Compression identifiers:

      <table>
        <tr> <th BGCOLOR="#EBEBEB">Column</th> <th BGCOLOR="#EBEBEB">Value</th> <th BGCOLOR="#EBEBEB">Description</th> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>5</td> <td>MSNumpress (linear) with lossy compression</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>6</td> <td>MSNumpress (short logged float) with lossy compression</td> </tr>
      </table>
  */
  class OPENMS_DLLAPI MSChromatogramParquetConsumer : public Interfaces::IMSDataConsumer
  {
  public:
    /**
      @brief Construct a parquet consumer for chromatogram export.

      @param[in] filename Output parquet filename.
      @param[in] run_id Run identifier to store with each chromatogram.
      @param[in] source_file Source mzML filename to store with each chromatogram.
      @param[in] transition_exp Transition metadata used to annotate chromatograms.
    */
    MSChromatogramParquetConsumer(const String& filename,
                                  UInt64 run_id,
                                  const String& source_file,
                                  const OpenSwath::LightTargetedExperiment& transition_exp);

    /// @brief Destructor flushes pending data and closes the parquet writer.
    ~MSChromatogramParquetConsumer() override;

    /// @brief Consume a spectrum (no-op; spectra are ignored for chromatogram export).
    void consumeSpectrum(SpectrumType& s) override;
    /// @brief Consume a chromatogram and append it to the parquet output.
    void consumeChromatogram(ChromatogramType& c) override;
    /**
      @brief Finalize and write the parquet file.

      Call this explicitly to surface write errors during normal control flow.
    */
    void finalize();
    /**
      @brief Reserve storage for expected data sizes.

      @param[in] expectedSpectra Expected number of spectra (ignored).
      @param[in] expectedChromatograms Expected number of chromatograms.
    */
    void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) override;
    /**
      @brief Set experimental settings (currently unused).

      @param[in] exp Experimental settings to store for context.
    */
    void setExperimentalSettings(const ExperimentalSettings& exp) override;

  private:
    std::unique_ptr<MSChromatogramParquetConsumerImpl> impl_;
  };
} // namespace OpenMS
