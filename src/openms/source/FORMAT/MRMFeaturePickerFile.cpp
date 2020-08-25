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

#include <OpenMS/FORMAT/MRMFeaturePickerFile.h>
#include <boost/regex.hpp>
#include <iostream>

namespace OpenMS
{
  void MRMFeaturePickerFile::load(
    const String& filename,
    std::vector<MRMFeaturePicker::ComponentParams>& cp_list,
    std::vector<MRMFeaturePicker::ComponentGroupParams>& cgp_list
  )
  {
    cp_list.clear();
    cgp_list.clear();
    CsvFile::load(filename, ',', false);
    StringList sl;
    std::map<String, Size> headers;
    if (rowCount() >= 2) // no need to read headers if that's the only line inside the file
    {
      getRow(0, sl);
      for (Size i = 0; i < sl.size(); ++i)
      {
        headers[sl[i]] = i; // for each header found, assign an index value to it
      }
      if (!(headers.count("component_name") && headers.count("component_group_name")))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Columns component_name and/or component_group_name not found.");
      }
    }
    for (Size i = 1; i < rowCount(); ++i)
    {
      getRow(i, sl);
      MRMFeaturePicker::ComponentParams cp;
      MRMFeaturePicker::ComponentGroupParams cgp;
      if (extractParamsFromLine_(sl, headers, cp, cgp))
      {
        cp_list.push_back(cp);
        // The following lines check if cgp is already present in cgp_list
        // If that is not the case, the extracted cgp is pushed
        bool cgp_found = std::any_of(
          cgp_list.begin(),
          cgp_list.end(),
          [&cgp](MRMFeaturePicker::ComponentGroupParams current_cgp)
          {
            return cgp.component_group_name == current_cgp.component_group_name;
          }
        );
        if (!cgp_found)
        {
          cgp_list.push_back(cgp);
        }
      }
    }
  }

  bool MRMFeaturePickerFile::extractParamsFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    MRMFeaturePicker::ComponentParams& cp,
    MRMFeaturePicker::ComponentGroupParams& cgp
  ) const
  {
    cp.component_name = line[headers.find("component_name")->second]; // save the component_name value
    cp.component_group_name = line[headers.find("component_group_name")->second]; // save the component_group_name value
    if (cp.component_name.empty() || cp.component_group_name.empty()) // component_name and component_group_name must not be empty
    {
      return false;
    }
    cgp.component_group_name = cp.component_group_name; // save the component_group_name also into cgp
    for (const std::pair<const String, Size>& h : headers) // parse the parameters
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("TransitionGroupPicker:(?!PeakPickerMRM:)(.+)")))
      {
        setCastValue_(String(m[1]), line[i], cgp.params);
      }
      else if (boost::regex_search(header, m, boost::regex("TransitionGroupPicker:PeakPickerMRM:(.+)")))
      {
        setCastValue_(String(m[1]), line[i], cp.params);
      }
    }
    return true;
  }

  void MRMFeaturePickerFile::setCastValue_(const String& key, const String& value, Param& params) const
  {
    if (value.empty()) // if the value is empty, don't set it
    {
      return;
    }
    const std::vector<String> double_headers = {
      "gauss_width", "peak_width", "signal_to_noise", "sn_win_len", "stop_after_intensity_ratio",
      "min_peak_width", "recalculate_peaks_max_z", "minimal_quality", "resample_boundary"
    };
    const std::vector<String> bool_headers = {
      "use_gauss", "write_sn_log_messages", "remove_overlapping_peaks", "recalculate_peaks",
      "use_precursors", "compute_peak_quality", "compute_peak_shape_metrics"
    };
    const std::vector<String> uint_headers = {
      "sgolay_frame_length", "sgolay_polynomial_order", "sn_bin_count"
    };
    const std::vector<String> int_headers = {
      "stop_after_feature"
    };
    if (std::find(double_headers.begin(), double_headers.end(), key) != double_headers.end())
    {
      params.setValue(key, value.toDouble());
    }
    else if (std::find(bool_headers.begin(), bool_headers.end(), key) != bool_headers.end())
    {
      params.setValue(key, value == "true" || value == "TRUE" ? "true" : "false");
    }
    else if (std::find(uint_headers.begin(), uint_headers.end(), key) != uint_headers.end())
    {
      params.setValue(key, static_cast<UInt>(value.toDouble()));
    }
    else if (std::find(int_headers.begin(), int_headers.end(), key) != int_headers.end())
    {
      params.setValue(key, value.toInt());
    }
    else // no conversion for class' parameters of type String
    {
      params.setValue(key, value);
    }
  }
}
