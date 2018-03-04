// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <boost/regex.hpp>

namespace OpenMS
{
  void MRMFeatureQCFile::load(const String& filename, MRMFeatureQC& mrmfqc) const
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
    if (headers.count("component_name")) // load component file
    {
      mrmfqc.component_qcs.clear();
      for (Size i = 1; i < csv.rowCount(); ++i)
      {
        csv.getRow(i, sl);
        pushValuesFromLine_(sl, headers, mrmfqc.component_qcs);
      }
    }
    else if (headers.count("component_group_name")) // load component group file
    {
      mrmfqc.component_group_qcs.clear();
      for (Size i = 1; i < csv.rowCount(); ++i)
      {
        csv.getRow(i, sl);
        pushValuesFromLine_(sl, headers, mrmfqc.component_group_qcs);
      }
    }
    else
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The file must contain one of the two following columns: component_name, component_group_name.");
    }
  }

  void MRMFeatureQCFile::pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
  ) const
  {
    MRMFeatureQC::ComponentQCs c;
    std::map<String, Size>::const_iterator it;
    it = headers.find("component_name");
    c.component_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("retention_time_l");
    c.retention_time_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("retention_time_u");
    c.retention_time_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("intensity_l");
    c.intensity_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("intensity_u");
    c.intensity_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("overall_quality_l");
    c.overall_quality_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("overall_quality_u");
    c.overall_quality_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    for (const std::pair<String, Size>& h : headers) // parse the parameters
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("metaValue_(.+)_(l|u)")))
      {
        setPairValue_(String(m[1]), line[i], String(m[2]), c.meta_value_qc);
      }
    }
    c_qcs.push_back(c); // TODO: check for any condition before pushing?
  }

  void MRMFeatureQCFile::pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
  ) const
  {
    MRMFeatureQC::ComponentGroupQCs cg;
    std::map<String, Size>::const_iterator it;
    it = headers.find("component_group_name");
    cg.component_group_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("retention_time_l");
    cg.retention_time_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("retention_time_u");
    cg.retention_time_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("intensity_l");
    cg.intensity_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("intensity_u");
    cg.intensity_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("overall_quality_l");
    cg.overall_quality_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("overall_quality_u");
    cg.overall_quality_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("n_heavy_l");
    cg.n_heavy_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_heavy_u");
    cg.n_heavy_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_light_l");
    cg.n_light_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_light_u");
    cg.n_light_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_detecting_l");
    cg.n_detecting_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_detecting_u");
    cg.n_detecting_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_quantifying_l");
    cg.n_quantifying_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_quantifying_u");
    cg.n_quantifying_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_identifying_l");
    cg.n_identifying_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_identifying_u");
    cg.n_identifying_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_transitions_l");
    cg.n_transitions_l = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("n_transitions_u");
    cg.n_transitions_u = it != headers.end() ? line[it->second].toInt() : 0;
    it = headers.find("ion_ratio_pair_name_1");
    cg.ion_ratio_pair_name_1 = it != headers.end() ? line[it->second] : "";
    it = headers.find("ion_ratio_pair_name_2");
    cg.ion_ratio_pair_name_2 = it != headers.end() ? line[it->second] : "";
    it = headers.find("ion_ratio_l");
    cg.ion_ratio_l = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("ion_ratio_u");
    cg.ion_ratio_u = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("ion_ratio_feature_name");
    cg.ion_ratio_feature_name = it != headers.end() ? line[it->second] : "";

    for (const std::pair<String, Size>& h : headers) // parse the parameters
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("metaValue_(.+)_(l|u)")))
      {
        setPairValue_(String(m[1]), line[i], String(m[2]), cg.meta_value_qc);
      }
    }
    cg_qcs.push_back(cg); // TODO: check for any condition before pushing?
  }

  void MRMFeatureQCFile::setPairValue_(
    const String& key,
    const String& value,
    const String& boundary,
    std::map<String, std::pair<double,double>>& meta_values_qc
  ) const
  {
    std::map<String, std::pair<double,double>>::iterator it = meta_values_qc.find(key);
    if (it != meta_values_qc.end())
    {
      if (boundary == "l") it->second.first = value.toDouble();
      else it->second.second = value.toDouble();
    }
    else
    {
      meta_values_qc[key] = boundary == "l"
        ? std::make_pair(value.toDouble(), 0.0)
        : std::make_pair(0.0, value.toDouble());
    }
  }

  //TODO
  // void MRMFeatureQCFile::store(const String & filename, const MRMFeatureQC & mrmfqc)
  // {
  //   // TODO: pending fix to CsvFile::fstore()
  // }

} // namespace OpenMS
