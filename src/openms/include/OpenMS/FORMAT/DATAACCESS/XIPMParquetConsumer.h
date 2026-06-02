// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <memory>

namespace OpenMS
{
  class XIPMParquetConsumerImpl;

  /**
    @brief Writes targeted mz/RT/IM peak maps to a Parquet file (.xipm).

    The schema includes run metadata, precursor/transition annotations,
    extraction window metadata, compressed mz/RT/ion mobility/intensity arrays,
    and compression flags.

    The Parquet output has the following columns (one row per extracted peak
    map):

      <table>
        <tr> <th BGCOLOR="#EBEBEB">Column</th> <th BGCOLOR="#EBEBEB">Type</th> <th BGCOLOR="#EBEBEB">Description</th> </tr>
        <tr> <td>RUN_ID</td> <td>int64</td> <td>Run identifier</td> </tr>
        <tr> <td>SOURCE_FILE</td> <td>string (nullable)</td> <td>Input source filename</td> </tr>
        <tr> <td>MS_LEVEL</td> <td>int64</td> <td>MS level (1 for precursor peak maps, 2 for fragment peak maps)</td> </tr>
        <tr> <td>PEAKMAP_TYPE</td> <td>string</td> <td>Peak-map type (<tt>precursor</tt> or <tt>transition</tt>)</td> </tr>
        <tr> <td>PRECURSOR_ID</td> <td>int64 (nullable)</td> <td>Precursor id</td> </tr>
        <tr> <td>TRANSITION_ID</td> <td>int64 (nullable)</td> <td>Transition id</td> </tr>
        <tr> <td>MODIFIED_SEQUENCE</td> <td>string (nullable)</td> <td>Modified peptide sequence</td> </tr>
        <tr> <td>PRECURSOR_CHARGE</td> <td>int64 (nullable)</td> <td>Precursor charge</td> </tr>
        <tr> <td>PRODUCT_CHARGE</td> <td>int64 (nullable)</td> <td>Product charge</td> </tr>
        <tr> <td>DETECTING_TRANSITION</td> <td>int64 (nullable)</td> <td>Detecting transition flag</td> </tr>
        <tr> <td>PRECURSOR_DECOY</td> <td>int64 (nullable)</td> <td>Precursor decoy flag</td> </tr>
        <tr> <td>PRODUCT_DECOY</td> <td>int64 (nullable)</td> <td>Product decoy flag</td> </tr>
        <tr> <td>TRANSITION_ORDINAL</td> <td>int64 (nullable)</td> <td>Transition ordinal</td> </tr>
        <tr> <td>TRANSITION_TYPE</td> <td>string (nullable)</td> <td>Transition type (e.g. y, b)</td> </tr>
        <tr> <td>ANNOTATION</td> <td>string (nullable)</td> <td>Transition annotation (e.g. y3^1) or precursor annotation</td> </tr>
        <tr> <td>TARGET_MZ</td> <td>float64</td> <td>Extraction target m/z</td> </tr>
        <tr> <td>TARGET_RT</td> <td>float64 (nullable)</td> <td>Extraction target RT in the run coordinate system</td> </tr>
        <tr> <td>TARGET_ION_MOBILITY</td> <td>float64 (nullable)</td> <td>Extraction target ion mobility</td> </tr>
        <tr> <td>RT_START</td> <td>float64 (nullable)</td> <td>Extraction window start in RT</td> </tr>
        <tr> <td>RT_END</td> <td>float64 (nullable)</td> <td>Extraction window end in RT</td> </tr>
        <tr> <td>MZ_DATA</td> <td>binary</td> <td>Compressed mz array</td> </tr>
        <tr> <td>RT_DATA</td> <td>binary</td> <td>Compressed RT array</td> </tr>
        <tr> <td>MOBILITY_DATA</td> <td>binary</td> <td>Compressed ion mobility array</td> </tr>
        <tr> <td>INTENSITY_DATA</td> <td>binary</td> <td>Compressed intensity array</td> </tr>
        <tr> <td>MZ_COMPRESSION</td> <td>int64</td> <td>m/z compression scheme id</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>int64</td> <td>RT compression scheme id</td> </tr>
        <tr> <td>MOBILITY_COMPRESSION</td> <td>int64</td> <td>Ion mobility compression scheme id</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>int64</td> <td>Intensity compression scheme id</td> </tr>
      </table>

    Compression identifiers:

      <table>
        <tr> <th BGCOLOR="#EBEBEB">Column</th> <th BGCOLOR="#EBEBEB">Value</th> <th BGCOLOR="#EBEBEB">Description</th> </tr>
        <tr> <td>MZ_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>MZ_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>MZ_COMPRESSION</td> <td>5</td> <td>MSNumpress (linear) with lossy compression</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>RT_COMPRESSION</td> <td>5</td> <td>MSNumpress (linear) with lossy compression</td> </tr>
        <tr> <td>MOBILITY_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>MOBILITY_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>MOBILITY_COMPRESSION</td> <td>5</td> <td>MSNumpress (linear) with lossy compression</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>0</td> <td>No compression (raw doubles)</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>1</td> <td>Zlib-compressed raw doubles</td> </tr>
        <tr> <td>INTENSITY_COMPRESSION</td> <td>6</td> <td>MSNumpress (short logged float) with lossy compression</td> </tr>
      </table>
  */
  class OPENMS_DLLAPI XIPMParquetConsumer
  {
  public:
    /**
      @brief Construct a parquet writer for targeted peak-map export.

      @param[in] filename Output parquet filename
      @param[in] transition_exp Transition metadata used to annotate rows
    */
    XIPMParquetConsumer(const String& filename,
                        const OpenSwath::LightTargetedExperiment& transition_exp);

    /// Destructor flushes pending rows and closes the parquet writer.
    ~XIPMParquetConsumer();

    /**
      @brief Consume one targeted peak map.

      @param[in] peak_map Extracted peak map to append
      @param[in] run_id Run identifier
      @param[in] source_file Source mzML / raw file name
      @param[in] ms_level MS level of the extraction (1 = precursor, 2 = fragment)
    */
    void consumePeakMap(const PeakMapExtractor::ExtractedPeakMap& peak_map,
                        UInt64 run_id,
                        const String& source_file,
                        Int64 ms_level);

    /// Finalize and write the parquet file.
    void finalize();

    /**
      @brief Reserve storage for the expected number of peak maps.

      @param[in] expectedPeakMaps The expected number of rows
    */
    void setExpectedSize(Size expectedPeakMaps);

  private:
    std::unique_ptr<XIPMParquetConsumerImpl> impl_;
  };
} // namespace OpenMS
