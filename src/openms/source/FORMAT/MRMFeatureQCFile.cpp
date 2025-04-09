// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <boost/regex.hpp>

namespace OpenMS
{
  void MRMFeatureQCFile::load(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) const
  {
    CsvFile csv(filename, ',', false, -1);
    StringList sl;
    std::map<String, Size> headers;
    if (csv.rowCount() > 0) // avoid accessing a row in an empty file
    {
      csv.getRow(0, sl);
    }
    for (Size i = 0; i < sl.size(); ++i)
    {
      headers[sl[i]] = i; // for each header found, assign an index value to it
    }
    if (!is_component_group) // load component file
    {
      mrmfqc.component_qcs.clear();
      for (Size i = 1; i < csv.rowCount(); ++i)
      {
        csv.getRow(i, sl);
        pushValuesFromLine_(sl, headers, mrmfqc.component_qcs);
      }
    }
    else // load component group file
    {
      mrmfqc.component_group_qcs.clear();
      for (Size i = 1; i < csv.rowCount(); ++i)
      {
        csv.getRow(i, sl);
        pushValuesFromLine_(sl, headers, mrmfqc.component_group_qcs);
      }
    }
  }

  void MRMFeatureQCFile::pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
  ) const
  {
    MRMFeatureQC::ComponentQCs c;
    c.component_name = getCastValue_(headers, line, "component_name", "");
    if (c.component_name.empty()) return;
    c.retention_time_l = getCastValue_(headers, line, "retention_time_l", 0.0);
    c.retention_time_u = getCastValue_(headers, line, "retention_time_u", 1e12);
    c.intensity_l = getCastValue_(headers, line, "intensity_l", 0.0);
    c.intensity_u = getCastValue_(headers, line, "intensity_u", 1e12);
    c.overall_quality_l = getCastValue_(headers, line, "overall_quality_l", 0.0);
    c.overall_quality_u = getCastValue_(headers, line, "overall_quality_u", 1e12);
    for (const std::pair<const String, Size>& h : headers) // parse the parameters
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("metaValue_(.+)_(l|u)"))) // capture the metavalue name and the boundary and save them to m[1] and m[2]
      {
        setPairValue_(String(m[1]), line[i], String(m[2]), c.meta_value_qc);
      }
    }
    c_qcs.push_back(c);
  }

  void MRMFeatureQCFile::pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
  ) const
  {
    MRMFeatureQC::ComponentGroupQCs cg;
    cg.component_group_name = getCastValue_(headers, line, "component_group_name", "");
    if (cg.component_group_name.empty())
    {
      return;
    }
    cg.retention_time_l = getCastValue_(headers, line, "retention_time_l", 0.0);
    cg.retention_time_u = getCastValue_(headers, line, "retention_time_u", 1e12);
    cg.intensity_l = getCastValue_(headers, line, "intensity_l", 0.0);
    cg.intensity_u = getCastValue_(headers, line, "intensity_u", 1e12);
    cg.overall_quality_l = getCastValue_(headers, line, "overall_quality_l", 0.0);
    cg.overall_quality_u = getCastValue_(headers, line, "overall_quality_u", 1e12);
    cg.n_heavy_l = getCastValue_(headers, line, "n_heavy_l", 0);
    cg.n_heavy_u = getCastValue_(headers, line, "n_heavy_u", 100);
    cg.n_light_l = getCastValue_(headers, line, "n_light_l", 0);
    cg.n_light_u = getCastValue_(headers, line, "n_light_u", 100);
    cg.n_detecting_l = getCastValue_(headers, line, "n_detecting_l", 0);
    cg.n_detecting_u = getCastValue_(headers, line, "n_detecting_u", 100);
    cg.n_quantifying_l = getCastValue_(headers, line, "n_quantifying_l", 0);
    cg.n_quantifying_u = getCastValue_(headers, line, "n_quantifying_u", 100);
    cg.n_identifying_l = getCastValue_(headers, line, "n_identifying_l", 0);
    cg.n_identifying_u = getCastValue_(headers, line, "n_identifying_u", 100);
    cg.n_transitions_l = getCastValue_(headers, line, "n_transitions_l", 0);
    cg.n_transitions_u = getCastValue_(headers, line, "n_transitions_u", 100);
    cg.ion_ratio_pair_name_1 = getCastValue_(headers, line, "ion_ratio_pair_name_1", "");
    cg.ion_ratio_pair_name_2 = getCastValue_(headers, line, "ion_ratio_pair_name_2", "");
    cg.ion_ratio_l = getCastValue_(headers, line, "ion_ratio_l", 0.0);
    cg.ion_ratio_u = getCastValue_(headers, line, "ion_ratio_u", 1e12);
    cg.ion_ratio_feature_name = getCastValue_(headers, line, "ion_ratio_feature_name", "");
    for (const std::pair<const String, Size>& h : headers) // parse the parameters
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("metaValue_(.+)_(l|u)"))) // capture the metavalue name and the boundary and save them to m[1] and m[2]
      {
        setPairValue_(String(m[1]), line[i], String(m[2]), cg.meta_value_qc);
      }
    }
    cg_qcs.push_back(cg);
  }

  void MRMFeatureQCFile::setPairValue_(
    const String& key,
    const String& value,
    const String& boundary,
    std::map<String, std::pair<double,double>>& meta_values_qc
  ) const
  {
    std::map<String, std::pair<double,double>>::iterator it = meta_values_qc.find(key);
    const double cast_value = value.empty() ? (boundary == "l" ? 0.0 : 1e12) : std::stod(value);
    if (it != meta_values_qc.end())
    {
      if (boundary == "l") it->second.first = cast_value;
      else it->second.second = cast_value;
    }
    else
    {
      meta_values_qc[key] = boundary == "l"
        ? std::make_pair(cast_value, 1e12)
        : std::make_pair(0.0, cast_value);
    }
  }

  Int MRMFeatureQCFile::getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const Int default_value
  ) const
  {
    std::map<String, Size>::const_iterator it = headers.find(header);
    return it != headers.end() && !line[it->second].empty()
      ? std::stoi(line[it->second])
      : default_value;
  }

  double MRMFeatureQCFile::getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const double default_value
  ) const
  {
    std::map<String, Size>::const_iterator it = headers.find(header);
    return it != headers.end() && !line[it->second].empty()
      ? std::stod(line[it->second])
      : default_value;
  }

  String MRMFeatureQCFile::getCastValue_(
    const std::map<String, Size>& headers,
    const StringList& line,
    const String& header,
    const String& default_value
  ) const
  {
    std::map<String, Size>::const_iterator it = headers.find(header);
    return it != headers.end() && !line[it->second].empty()
      ? line[it->second]
      : default_value;
  }

  void MRMFeatureQCFile::store(const String & filename, const MRMFeatureQC & mrmfqc, const bool is_component_group)
  {
    if (is_component_group) 
    {
      if (mrmfqc.component_group_qcs.empty()) return;

      // Store the ComponentGroupQCs
      clear(); // clear the buffer_

      // Make the ComponentGroupQCs headers
      StringList headers = { "component_group_name", "retention_time_l", "retention_time_u", "intensity_l", "intensity_u", "overall_quality_l", "overall_quality_u",
        "n_heavy_l", "n_heavy_u", "n_light_l", "n_light_u", "n_detecting_l", "n_detecting_u", "n_quantifying_l", "n_quantifying_u", "n_identifying_l", "n_identifying_u", "n_transitions_l", "n_transitions_u",
        "ion_ratio_pair_name_1", "ion_ratio_pair_name_2", "ion_ratio_l", "ion_ratio_u", "ion_ratio_feature_name" };
      for (const auto& meta_data : mrmfqc.component_group_qcs.at(0).meta_value_qc) 
      {
        headers.push_back("metaValue_" + meta_data.first + "_l");
        headers.push_back("metaValue_" + meta_data.first + "_u");
      }
      addRow(headers);

      // Make the ComponentGroupQCs rows
      for (const auto& component_qc : mrmfqc.component_group_qcs)
      {
        StringList row(headers.size());
        row[0] = component_qc.component_group_name;
        row[1] = component_qc.retention_time_l;
        row[2] = component_qc.retention_time_u;
        row[3] = component_qc.intensity_l;
        row[4] = component_qc.intensity_u;
        row[5] = component_qc.overall_quality_l;
        row[6] = component_qc.overall_quality_u;
        row[7] = component_qc.n_heavy_l;
        row[8] = component_qc.n_heavy_u;
        row[9] = component_qc.n_light_l;
        row[10] = component_qc.n_light_u;
        row[11] = component_qc.n_detecting_l;
        row[12] = component_qc.n_detecting_u;
        row[13] = component_qc.n_quantifying_l;
        row[14] = component_qc.n_quantifying_u;
        row[15] = component_qc.n_identifying_l;
        row[16] = component_qc.n_identifying_u;
        row[17] = component_qc.n_transitions_l;
        row[18] = component_qc.n_transitions_u;
        row[19] = component_qc.ion_ratio_pair_name_1;
        row[20] = component_qc.ion_ratio_pair_name_2;
        row[21] = component_qc.ion_ratio_l;
        row[22] = component_qc.ion_ratio_u;
        row[23] = component_qc.ion_ratio_feature_name;
        size_t meta_data_iter = 24;
        for (const auto& meta_data : component_qc.meta_value_qc) {
          row[meta_data_iter] = meta_data.second.first;
          ++meta_data_iter;
          row[meta_data_iter] = meta_data.second.second;
          ++meta_data_iter;
        }
        addRow(row);
      }

      CsvFile::store(filename);
    } 
    else 
    {
      if (mrmfqc.component_qcs.empty()) return;

      // Store the ComponentQCs
      clear(); // clear the buffer_

      // Make the ComponentQCs headers
      StringList headers = { "component_name","retention_time_l","retention_time_u","intensity_l","intensity_u","overall_quality_l","overall_quality_u" };
      for (const auto& meta_data : mrmfqc.component_qcs.at(0).meta_value_qc) 
      {
        headers.push_back("metaValue_" + meta_data.first + "_l");
        headers.push_back("metaValue_" + meta_data.first + "_u");
      }
      addRow(headers);

      // Make the ComponentQCs rows
      for (const auto& component_qc : mrmfqc.component_qcs)
      {
        StringList row(headers.size());
        row[0] = component_qc.component_name;
        row[1] = component_qc.retention_time_l;
        row[2] = component_qc.retention_time_u;
        row[3] = component_qc.intensity_l;
        row[4] = component_qc.intensity_u;
        row[5] = component_qc.overall_quality_l;
        row[6] = component_qc.overall_quality_u;
        size_t meta_data_iter = 7;
        for (const auto& meta_data : component_qc.meta_value_qc) 
        {
          row[meta_data_iter] = meta_data.second.first;
          ++meta_data_iter;
          row[meta_data_iter] = meta_data.second.second;
          ++meta_data_iter;
        }
        addRow(row);
      }
      CsvFile::store(filename);

    }
  }

} // namespace OpenMS
