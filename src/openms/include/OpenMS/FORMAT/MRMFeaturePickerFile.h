// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
  /**
    @brief _MRMFeaturePickerFile_ loads components and components groups parameters
    from a .csv file.

    The structures defined in [MRMFeaturePicker](@ref MRMFeaturePicker) are used.

    It is required that columns `component_name` and `component_group_name` are present.
    Lines whose `component_name`'s or `component_group_name`'s value is an empty string, will be skipped.
    The class supports the absence of information within other columns.

    A reduced example of the expected format (fewer columns are shown here):
    > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
    > arg-L.arg-L_1.Heavy,arg-L,2,15
    > arg-L.arg-L_1.Light,arg-L,2,17
    > orn.orn_1.Heavy,orn,3,21
    > orn.orn_1.Light,orn,3,13
  */
  class OPENMS_DLLAPI MRMFeaturePickerFile :
    public CsvFile
  {
public:
    /// Constructor
    MRMFeaturePickerFile() = default;
    /// Destructor
    ~MRMFeaturePickerFile() override = default;

    /**
      @brief Loads the file's data and saves it into vectors of `ComponentParams` and `ComponentGroupParams`.

      The file is expected to contain at least two columns: `component_name` and `component_group_name`. Otherwise,
      an exception is thrown.

      If a component group (identified by its name) is found multiple times, only the first one is saved.

      @param[in] filename Path to the .csv input file
      @param[out] cp_list Component params are saved in this list
      @param[out] cgp_list Component Group params are saved in this list

      @throw Exception::MissingInformation If the required columns are not found.
      @throw Exception::FileNotFound If input file is not found.
    */
    void load(
      const String& filename,
      std::vector<MRMFeaturePicker::ComponentParams>& cp_list,
      std::vector<MRMFeaturePicker::ComponentGroupParams>& cgp_list
    );

protected:
    /**
      @brief Extracts the information from a `StringList` and saves it into the correct data structures.

      @param[in] line The line parsed from the input file
      @param[in] headers A mapping from a given header to its value's position
      @param[out] cp The extracted component parameters
      @param[out] cgp The extracted component group parameters

      @return Returns `false` if `component_name` or `component_group_name` are empty strings. Otherwise, it returns `true`.
    */
    bool extractParamsFromLine_(
      const StringList& line,
      const std::map<String, Size>& headers,
      MRMFeaturePicker::ComponentParams& cp,
      MRMFeaturePicker::ComponentGroupParams& cgp
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
}

