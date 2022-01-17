// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
     @param map MzTabMMetaData
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
     @param map MzTabMMetaData
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
     @param map MzTabMMetaData
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
