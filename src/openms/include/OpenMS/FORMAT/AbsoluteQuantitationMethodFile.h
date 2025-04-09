// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
  /**
    @brief File adapter for AbsoluteQuantitationMethod files.

    Loads and stores .csv or .tsv files describing an AbsoluteQuantitationMethod.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI AbsoluteQuantitationMethodFile :
    private CsvFile
  {
public:
    ///Default constructor
    AbsoluteQuantitationMethodFile() = default;
    ///Destructor
    ~AbsoluteQuantitationMethodFile() override = default;

    /**
      @brief Loads an AbsoluteQuantitationMethod file.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing

      @param[in] filename The input file name.
      @param[out] aqm_list Output variable where the AbsoluteQuantitationMethod data is loaded.
    */
    void load(const String & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list);

    /**
      @brief Stores an AbsoluteQuantitationMethod file.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created

      @param[in] filename The output file name.
      @param[in] aqm_list The AbsoluteQuantitationMethod data to write into the file.
    */
    void store(const String & filename, const std::vector<AbsoluteQuantitationMethod> & aqm_list);

protected:
    /**
      @brief Parses a line into the members of AbsoluteQuantitationMethod.

      @param[in] line A line of the .csv file.
      @param[in] headers A map of header strings to column positions.
      @param[out] aqm AbsoluteQuantitationMethod.
    */
    void parseLine_(
      const StringList & line,
      const std::map<String, Size> & headers,
      AbsoluteQuantitationMethod & aqm
    ) const;

    /**
      @brief Helper method which takes care of converting the given value to the desired type,
      based on the header (here `key`) information.

      @param[in] key The header name with which the correct conversion is chosen
      @param[in] value The value to be converted
      @param[in,out] params The object where the new value is saved
    */
    void setCastValue_(const String& key, const String& value, Param& params) const;

  };

} // namespace OpenMS

