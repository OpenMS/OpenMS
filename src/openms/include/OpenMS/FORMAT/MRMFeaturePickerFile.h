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
    > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length
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
    ~MRMFeaturePickerFile() = default;

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

