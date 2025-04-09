// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

namespace OpenMS
{
  /**
      @brief File adapter for MRMFeatureQC files.

      Loads and stores .csv or .tsv files describing an MRMFeatureQC.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MRMFeatureQCFile :
    private CsvFile,
    public ProgressLogger
  {
public:
  /// Default constructor
  MRMFeatureQCFile() = default;
  /// Destructor
  ~MRMFeatureQCFile() override = default;

  /**
    @brief Loads an MRMFeatureQC file.

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing

    @param[in] filename The path to the input file
    @param[in,out] mrmfqc The output class which will contain the criteria
    @param[in] is_component_group true if the user intends to load ComponentGroupQCs data, false otherwise
  */
  void load(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) const;

  /*
    @brief Stores an MRMFeatureQC file.

    @exception Exception::UnableToCreateFile is thrown if the file could not be created

    @param[in] filename The path to the input file
    @param[in] mrmfqc The output class which will contain the criteria
    @param[in] is_component_group true if the user intends to store ComponentGroupQCs data, false otherwise
  */
  void store(const String& filename, const MRMFeatureQC& mrmfqc, const bool is_component_group);
protected:
  /**
    @brief Save values from a line to a `ComponentQCs`.

    @note Lines missing the `component_name` value will be skipped.

    @param[in] line A line containing the values from a row in the input file
    @param[in] headers A mapping from headers names to position indices
    @param[out] c_qcs The output will be saved in a new element of this vector
  */
  void pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
  ) const;

  /**
    @brief Save values from a line to a `ComponentGroupQCs`.

    @note Lines missing the `component_group_name` value will be skipped.

    @param[in] line A line containing the values from a row in the input file
    @param[in] headers A mapping from headers names to position indices
    @param[out] cg_qcs The output will be saved in a new element of this vector
  */
  void pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
  ) const;

  /**
    @brief Set one of the values in a pair

    The method is given in input a map from Strings to pairs.
    Assuming the key is present within the map, its mapped value will be updated
    with the given `value`, in the correct `boundary` position.
    If the key is not found, a new pair will be created and both its values set
    (in this case, a default value of 0.0 is given to the other pair's element).

    @param[in] key The metavalue name
    @param[in] value The mapped pair
    @param[in] boundary "l" for lower bound or "u" for upper bound
    @param[out] meta_values_qc The map containing the metavalues and pairs
  */
  void setPairValue_(
    const String& key,
    const String& value,
    const String& boundary,
    std::map<String, std::pair<double,double>>& meta_values_qc
  ) const;

  /**
    @brief Extracts a column's value from a line.

    The method looks for the value found within `line[headers[header]]`.
    If the information is present and its value is valid, it will be converted
    to `Int` and returned.
    Otherwise, `default_value` (provided by the user) is returned.

    @param[in] headers The mapping from columns' name to positions' indices
    @param[in] line A list of strings containing a single row's values
    @param[in] header The desired value's column name
    @param[in] default_value A default value to return in case the information is not found or invalid
    @return The found information (if found and valid) converted to `Int`. Otherwise `default_value`.
  */
  Int getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const Int default_value
  ) const;

  /**
    @brief Extracts a column's value from a line.

    The method looks for the value found within `line[headers[header]]`.
    If the information is present and its value is valid, it will be converted
    to `double` and returned.
    Otherwise, `default_value` (provided by the user) is returned.

    @param[in] headers The mapping from columns' name to positions' indices
    @param[in] line A list of strings containing a single row's values
    @param[in] header The desired value's column name
    @param[in] default_value A default value to return in case the information is not found or invalid
    @return The found information (if found and valid) converted to `double`. Otherwise `default_value`.
  */
  double getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const double default_value
  ) const;

  /**
    @brief Extracts a column's value from a line.

    The method looks for the value found within `line[headers[header]]`.
    If the information is present and its value is valid, it will be converted
    to `String` and returned.
    Otherwise, `default_value` (provided by the user) is returned.

    @param[in] headers The mapping from columns' name to positions' indices
    @param[in] line A list of strings containing a single row's values
    @param[in] header The desired value's column name
    @param[in] default_value A default value to return in case the information is not found or invalid
    @return The found information (if found and valid) converted to `String`. Otherwise `default_value`.
  */
  String getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const String& default_value
  ) const;

  };

} // namespace OpenMS

