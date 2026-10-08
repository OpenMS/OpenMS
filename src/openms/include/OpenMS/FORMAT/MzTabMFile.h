// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/MzTabM.h>

namespace OpenMS
{
  class SVOutStream;

  /**
    @brief File adapter for MzTab-M files

    @ingroup FileIO
  */

  class OPENMS_DLLAPI MzTabMFile
  {
  public:
    /// Default Constructor
    MzTabMFile();

    /// Default Destructor
    ~MzTabMFile();

    /// Store MzTabM file
    void store(const std::string& filename, const MzTabM& mztab_m) const;

  protected:

    /**
     @brief Generates the MzTabM MetaData Section
     @param[in] map MzTabMMetaData
     @param[out] sl Fill Stringlist with MztabM MetaData entries
    */
    void generateMzTabMMetaDataSection_(const MzTabMMetaData& map, StringList& sl) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param[in] meta MzTabMMetaData
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns in the header
     @return StringList with SMH entries
    */
    std::string generateMzTabMSmallMoleculeHeader_(const MzTabMMetaData& meta, const std::vector<std::string>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Section
     @param[in] row MzTabMSmallMoleculeSectionRow
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns per row
     @return StringList with SML entries
    */
    std::string generateMzTabMSmallMoleculeSectionRow_(const MzTabMSmallMoleculeSectionRow& row, const std::vector<std::string>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param[in] meta MzTabMMetaData
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns in the header
     @return StringList with SFH entries
    */
    std::string generateMzTabMSmallMoleculeFeatureHeader_(const MzTabMMetaData& meta, const std::vector<std::string>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Feature Section
     @param[in] row MzTabMSmallMoleculeFeatureSectionRow
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns per row
     @return StringList with SMF entries
    */
    std::string generateMzTabMSmallMoleculeFeatureSectionRow_(const MzTabMSmallMoleculeFeatureSectionRow& row, const std::vector<std::string>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param[in] meta MzTabMMetaData
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns in the header
     @return StringList with SEH entries
    */
    std::string generateMzTabMSmallMoleculeEvidenceHeader_(const MzTabMMetaData& meta, const std::vector<std::string>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Evidence Section
     @param[in] row MzTabMSmallMoleculeFeatureSectionRow
     @param[in] optional_columns Add optional columns
     @param[out] n_columns Stores the number of columns per row
     @return StringList with SME entries
    */
    std::string generateMzTabMSmallMoleculeEvidenceSectionRow_(const MzTabMSmallMoleculeEvidenceSectionRow& row, const std::vector<std::string>& optional_columns, size_t& n_columns) const;
  };

} // namespace OpenMS
