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
     * @param[in] filename The output file
     * @param[in] id_data The IdentificationData object
     */
    void store(const std::string& filename, const IdentificationData& id_data);

    /** @brief Write out a feature map to SQL-based OMS file
     *
     * @param[in] filename The output file
     * @param[in] features The feature map
     */
    void store(const std::string& filename, const FeatureMap& features);

    /** @brief Write out a consensus map to SQL-based OMS file
     *
     * @param[in] filename The output file
     * @param[in] consensus The consensus map
     */
    void store(const std::string& filename, const ConsensusMap& consensus);

    /** @brief Read in an OMS file and construct an IdentificationData object
     *
     * @param[out] filename The input file
     * @param[in] id_data The IdentificationData object
     */
    void load(const std::string& filename, IdentificationData& id_data);

    /** @brief Read in an OMS file and construct a feature map
     *
     * @param[out] filename The input file
     * @param[in] features The feature map
     */
    void load(const std::string& filename, FeatureMap& features);

    /** @brief Read in an OMS file and construct a consensus map
     *
     * @param[out] filename The input file
     * @param[in] consensus The consensus map
     */
    void load(const std::string& filename, ConsensusMap& consensus);

    /** @brief Read in an OMS file and write out the contents in JSON format
     *
     * @param[in] filename_in The input file (OMS)
     * @param[out] filename_out The output file (JSON)
     */
    void exportToJSON(const std::string& filename_in, const std::string& filename_out);

  protected:
    LogType log_type_;
  };
} // namespace OpenMS
