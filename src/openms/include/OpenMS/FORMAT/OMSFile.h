// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace OpenMS
{
  class FeatureMap;
  class ConsensusMap;

  /**
      @brief This class supports reading and writing of OMS files.

      OMS files are SQLite databases consisting of several tables.
  */
  class OPENMS_DLLAPI OMSFile: public ProgressLogger
  {
  public:
    /// Constructor (with option to set log type)
    explicit OMSFile(LogType log_type = LogType::NONE):
      log_type_(log_type)
    {
      setLogType(log_type);
    }

    /** @brief Write out an IdentificationData object to SQL-based OMS file
     *
     * @param filename The output file
     * @param id_data The IdentificationData object
     */
    void store(const String& filename, const IdentificationData& id_data);

    /** @brief Write out a feature map to SQL-based OMS file
     *
     * @param filename The output file
     * @param features The feature map
     */
    void store(const String& filename, const FeatureMap& features);

    /** @brief Write out a consensus map to SQL-based OMS file
     *
     * @param filename The output file
     * @param consensus The consensus map
     */
    void store(const String& filename, const ConsensusMap& consensus);

    /** @brief Read in an OMS file and construct an IdentificationData object
     *
     * @param filename The input file
     * @param id_data The IdentificationData object
     */
    void load(const String& filename, IdentificationData& id_data);

    /** @brief Read in an OMS file and construct a feature map
     *
     * @param filename The input file
     * @param features The feature map
     */
    void load(const String& filename, FeatureMap& features);

    /** @brief Read in an OMS file and construct a consensus map
     *
     * @param filename The input file
     * @param consensus The consensus map
     */
    void load(const String& filename, ConsensusMap& consensus);

    /** @brief Read in an OMS file and write out the contents in JSON format
     *
     * @param filename_in The input file (OMS)
     * @param filename_out The output file (JSON)
     */
    void exportToJSON(const String& filename_in, const String& filename_out);

  protected:
    LogType log_type_;
  };
} // namespace OpenMS
