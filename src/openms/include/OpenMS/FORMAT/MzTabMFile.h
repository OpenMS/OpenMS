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
  class String;
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
    void store(const String& filename, const MzTabM& mztab_m) const;

  protected:

    /**
     @brief Generates the MzTabM MetaData Section
     @param map MzTabMMetaData
     @param sl Fill Stringlist with MztabM MetaData entries
    */
    void generateMzTabMMetaDataSection_(const MzTabMMetaData& map, StringList& sl) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param meta MzTabMMetaData
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns in the header
     @return StringList with SMH entries
    */
    String generateMzTabMSmallMoleculeHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Section
     @param row MzTabMSmallMoleculeSectionRow
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns per row
     @return StringList with SML entries
    */
    String generateMzTabMSmallMoleculeSectionRow_(const MzTabMSmallMoleculeSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param meta MzTabMMetaData
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns in the header
     @return StringList with SFH entries
    */
    String generateMzTabMSmallMoleculeFeatureHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Feature Section
     @param row MzTabMSmallMoleculeFeatureSectionRow
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns per row
     @return StringList with SMF entries
    */
    String generateMzTabMSmallMoleculeFeatureSectionRow_(const MzTabMSmallMoleculeFeatureSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Header
     @param meta MzTabMMetaData
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns in the header
     @return StringList with SEH entries
    */
    String generateMzTabMSmallMoleculeEvidenceHeader_(const MzTabMMetaData& meta, const std::vector<String>& optional_columns, size_t& n_columns) const;

    /**
     @brief Generates the MzTabM Small Molecule Evidence Section
     @param row MzTabMSmallMoleculeFeatureSectionRow
     @param optional_columns Add optional columns
     @param n_columns Stores the number of columns per row
     @return StringList with SME entries
    */
    String generateMzTabMSmallMoleculeEvidenceSectionRow_(const MzTabMSmallMoleculeEvidenceSectionRow& row, const std::vector<String>& optional_columns, size_t& n_columns) const;
  };

} // namespace OpenMS
