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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MRMFeatureQCFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

namespace OpenMS
{

  MRMFeatureQCFile::MRMFeatureQCFile()
  {
  }

  MRMFeatureQCFile::~MRMFeatureQCFile()
  {
  }

  void MRMFeatureQCFile::load(const String & filename,
    MRMFeatureQC & mrmfqc)
  {
    // read in the .csv file
    char is = ',';
    bool ie = false; 
    Int first_n = -1;
    CsvFile::load(filename, is, ie, first_n);

    // parse the file
    std::map<String,int> headers;
    std::map<String,int> params_headers;
    StringList line, header;
    for (size_t i = 0; i < CsvFile::rowCount(); ++i)
    {
      if (i == 0) // header row
      {
        CsvFile::getRow(i, header);
        parseHeader_(header, headers, params_headers);
      }
      else
      {
        CsvFile::getRow(i, line);
        parseLine_(line, headers, params_headers, mrmfqc);
      }    
    }
  }

  void MRMFeatureQCFile::parseHeader_(StringList & line, std::map<String, int> & headers,
    std::map<String, int> & params_headers)
  {    
    // default header column positions
    headers["component_name"] = -1;
    headers["component_group_name"] = -1;
    headers["n_heavy_l"] = -1;
    headers["n_heavy_u"] = -1;
    headers["n_light_l"] = -1;
    headers["n_light_u"] = -1;
    headers["n_detecting_l"] = -1;
    headers["n_detecting_u"] = -1;
    headers["n_quantifying_l"] = -1;
    headers["n_quantifying_u"] = -1;
    headers["n_identifying_l"] = -1;
    headers["n_identifying_u"] = -1;
    headers["n_transitions_l"] = -1;
    headers["n_transitions_u"] = -1;
    headers["ion_ratio_pair_name_1"] = -1;
    headers["ion_ratio_pair_name_2"] = -1;
    headers["ion_ratio_l"] = -1;
    headers["ion_ratio_u"] = -1;
    headers["ion_ratio_feature_name"] = -1;
    headers["retention_time_l"] = -1;
    headers["retention_time_u"] = -1;
    headers["intensity_l"] = -1;
    headers["intensity_u"] = -1;
    headers["overall_quality_l"] = -1;
    headers["overall_quality_u"] = -1;
    String param_header = "metaValue_";
    
    // parse the header columns
    for (size_t i = 0; i < line.size(); ++i)
    {
      // parse transformation_model_params
      if (line[i].find(param_header) != String::npos) 
      {
        line[i].erase(0, param_header.size()); 
        params_headers[line[i]] = i;
      }      
      else // parse all other header entries
      {
        headers[line[i]] = i;
      }
    }
  }

  void MRMFeatureQCFile::parseLine_(StringList & line, std::map<String,int> & headers, 
    std::map<String,int> & params_headers,
    MRMFeatureQC & mrmfqc)
  {
    // component QCs
    MRMFeatureQC::ComponentQCs cqcs;
    cqcs.component_name = "";
    if (headers["component_name"] != -1)
    {
      cqcs.component_name = line[headers["component_name"]];
    }
    cqcs.retention_time_l = 0;
    if (headers["retention_time_l"] != -1)
    {
      cqcs.retention_time_l = (line[headers["retention_time_l"]].empty()) ? 0.0 : std::stod(line[headers["retention_time_l"]]);
    }
    cqcs.retention_time_u = 0;
    if (headers["retention_time_u"] != -1)
    {
      cqcs.retention_time_u = (line[headers["retention_time_u"]].empty()) ? 0.0 : std::stod(line[headers["retention_time_u"]]);
    }    
    cqcs.intensity_l = 0;
    if (headers["intensity_l"] != -1)
    {
      cqcs.intensity_l = (line[headers["intensity_l"]].empty()) ? 0.0 : std::stod(line[headers["intensity_l"]]);
    }
    cqcs.intensity_u = 0;
    if (headers["intensity_u"] != -1)
    {
      cqcs.intensity_u = (line[headers["intensity_u"]].empty()) ? 0.0 : std::stod(line[headers["intensity_u"]]);
    }
    cqcs.overall_quality_l = 0;
    if (headers["overall_quality_l"] != -1)
    {
      cqcs.overall_quality_l = (line[headers["overall_quality_l"]].empty()) ? 0.0 : std::stod(line[headers["overall_quality_l"]]);
    }
    cqcs.overall_quality_u = 0;
    if (headers["overall_quality_u"] != -1)
    {
      cqcs.overall_quality_u = (line[headers["overall_quality_u"]].empty()) ? 0.0 : std::stod(line[headers["overall_quality_u"]]);
    }
    // parse metaValues
    String meta_value_key = "";
    String lub = "";
    std::pair<double,double> lbub {0,0};
    for (auto const& kv : params_headers)
    {
      // split into meta_value_key and lub
      // example "meta_value_value_l" -> "meta_value_value" and "l"
      meta_value_key = kv.first.substr(0, kv.first.length()-2);
      lub = kv.first.substr(kv.first.length()-1, kv.first.length());
      if (cqcs.meta_value_qc.count(meta_value_key) == 0)
      {     
        cqcs.meta_value_qc[meta_value_key] = lbub;
      }
        
      // cast doubles
      if (lub == "l")
      {
        cqcs.meta_value_qc[meta_value_key].first = std::stod(line[kv.second]);
      }
      else if (lub == "u")
      {
        cqcs.meta_value_qc[meta_value_key].second = std::stod(line[kv.second]);
      }
      // cqcs.meta_value_qc
      
    }
    mrmfqc.component_qcs.push_back(cqcs);

    //component_group QCs
    MRMFeatureQC::ComponentGroupQCs cgqcs;
    cgqcs.component_group_name = "";
    if (headers["component_group_name"] != -1)
    {
      cgqcs.component_group_name = line[headers["component_group_name"]];
    }
    cgqcs.n_heavy_l = 0;
    if (headers["n_heavy_l"] != -1)
    {
      cgqcs.n_heavy_l = (line[headers["n_heavy_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_heavy_l"]]);
    }
    cgqcs.n_heavy_u = 0;
    if (headers["n_heavy_u"] != -1)
    {
      cgqcs.n_heavy_u = (line[headers["n_heavy_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_heavy_u"]]);
    }
    cgqcs.n_light_l = 0;
    if (headers["n_light_l"] != -1)
    {
      cgqcs.n_light_l = (line[headers["n_light_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_light_l"]]);
    }
    cgqcs.n_light_u = 0;
    if (headers["n_light_u"] != -1)
    {
      cgqcs.n_light_u = (line[headers["n_light_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_light_u"]]);
    } 
    cgqcs.n_detecting_l = 0;
    if (headers["n_detecting_l"] != -1)
    {
      cgqcs.n_detecting_l = (line[headers["n_detecting_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_detecting_l"]]);
    }
    cgqcs.n_detecting_u = 0;
    if (headers["n_detecting_u"] != -1)
    {
      cgqcs.n_detecting_u = (line[headers["n_detecting_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_detecting_u"]]);
    }
    cgqcs.n_quantifying_l = 0;
    if (headers["n_quantifying_l"] != -1)
    {
      cgqcs.n_quantifying_l = (line[headers["n_quantifying_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_quantifying_l"]]);
    }
    cgqcs.n_quantifying_u = 0;
    if (headers["n_quantifying_u"] != -1)
    {
      cgqcs.n_quantifying_u = (line[headers["n_quantifying_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_quantifying_u"]]);
    }
    cgqcs.n_identifying_l = 0;
    if (headers["n_identifying_l"] != -1)
    {
      cgqcs.n_identifying_l = (line[headers["n_identifying_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_identifying_l"]]);
    }
    cgqcs.n_identifying_u = 0;
    if (headers["n_identifying_u"] != -1)
    {
      cgqcs.n_identifying_u = (line[headers["n_identifying_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_identifying_u"]]);
    }
    cgqcs.n_transitions_l = 0;
    if (headers["n_transitions_l"] != -1)
    {
      cgqcs.n_transitions_l = (line[headers["n_transitions_l"]].empty()) ? 0.0 : std::stoi(line[headers["n_transitions_l"]]);
    }
    cgqcs.n_transitions_u = 0;
    if (headers["n_transitions_u"] != -1)
    {
      cgqcs.n_transitions_u = (line[headers["n_transitions_u"]].empty()) ? 0.0 : std::stoi(line[headers["n_transitions_u"]]);
    }
    cgqcs.ion_ratio_pair_name_1 = "";
    if (headers["ion_ratio_pair_name_1"] != -1)
    {
      cgqcs.ion_ratio_pair_name_1 = line[headers["ion_ratio_pair_name_1"]];
    }
    cgqcs.ion_ratio_pair_name_2 = "";
    if (headers["ion_ratio_pair_name_2"] != -1)
    {
      cgqcs.ion_ratio_pair_name_2 = line[headers["ion_ratio_pair_name_2"]];
    }
    cgqcs.ion_ratio_l = 0;
    if (headers["ion_ratio_l"] != -1)
    {
      cgqcs.ion_ratio_l = (line[headers["ion_ratio_l"]].empty()) ? 0.0 : std::stod(line[headers["ion_ratio_l"]]);
    }
    cgqcs.ion_ratio_u = 0;
    if (headers["ion_ratio_u"] != -1)
    {
      cgqcs.ion_ratio_u = (line[headers["ion_ratio_u"]].empty()) ? 0.0 : std::stod(line[headers["ion_ratio_u"]]);
    }
    cgqcs.ion_ratio_feature_name = "";
    if (headers["ion_ratio_feature_name"] != -1)
    {
      cgqcs.ion_ratio_feature_name = line[headers["ion_ratio_feature_name"]];
    }
    mrmfqc.component_group_qcs.push_back(cgqcs);
  }

  //TODO
  // void MRMFeatureQCFile::store(const String & filename, const MRMFeatureQC & mrmfqc)
  // {
  //   // TODO: pending fix to CsvFile::fstore()
  // }

} // namespace OpenMS
