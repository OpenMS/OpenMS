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

#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>

namespace OpenMS
{
  void AbsoluteQuantitationMethodFile::load(const String & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)
  {
    // read in the .csv file
    char is = ',';
    bool ie = false; 
    Int first_n = -1;
    CsvFile::load(filename, is, ie, first_n);

    // parse the file
    std::map<String, Int> headers;
    std::map<String, Int> params_headers;
    StringList line, header;
    AbsoluteQuantitationMethod aqm;
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
        parseLine_(line, headers, params_headers, aqm);
        aqm_list.push_back(aqm);
      }
    }
  }

  void AbsoluteQuantitationMethodFile::parseHeader_(
    StringList & line,
    std::map<String, Int> & headers,
    std::map<String, Int> & params_headers
  ) const
  {    
    // default header column positions
    headers["IS_name"] = -1;
    headers["component_name"] = -1;
    headers["feature_name"] = -1;
    headers["concentration_units"] = -1;
    headers["llod"] = -1;
    headers["ulod"] = -1;
    headers["lloq"] = -1;
    headers["uloq"] = -1;
    headers["correlation_coefficient"] = -1;
    headers["n_points"] = -1;
    headers["transformation_model"] = -1;
    String param_header = "transformation_model_param_";
    
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

  void AbsoluteQuantitationMethodFile::parseLine_(
    const StringList & line,
    const std::map<String, Int> & headers,
    const std::map<String, Int> & params_headers,
    AbsoluteQuantitationMethod & aqm
  ) const
  {
    // component, IS, and feature names
    String component_name = "";
    if (headers.at("component_name") != -1)
    {
      component_name = line[headers.at("component_name")];
    }
    aqm.setComponentName(component_name);
    String feature_name = "";
    if (headers.at("feature_name") != -1)
    {
      feature_name = line[headers.at("feature_name")];
    }
    aqm.setFeatureName(feature_name);
    String IS_name = "";
    if (headers.at("IS_name") != -1)
    {
      IS_name = line[headers.at("IS_name")];
    }
    aqm.setISName(IS_name);

    // LODs
    double llod = 0.0;
    if (headers.at("llod") != -1)
    {
      llod = (line[headers.at("llod")].empty()) ? 0.0 : std::stod(line[headers.at("llod")]);
    }
    aqm.setLLOD(llod);
    double ulod = 0.0;
    if (headers.at("ulod") != -1)
    {
      ulod = (line[headers.at("ulod")].empty()) ? 0.0 : std::stod(line[headers.at("ulod")]);
    }
    aqm.setULOD(ulod);

    // LOQs
    double lloq = 0.0;
    if (headers.at("lloq") != -1)
    {
      lloq = (line[headers.at("lloq")].empty()) ? 0.0 : std::stod(line[headers.at("lloq")]);
    }
    aqm.setLLOQ(lloq);
    double uloq = 0.0;
    if (headers.at("uloq") != -1)
    {
      uloq = (line[headers.at("uloq")].empty()) ? 0.0 : std::stod(line[headers.at("uloq")]);
    }
    aqm.setULOQ(uloq);

    // concentration units
    String concentration_units = "";
    if (headers.at("concentration_units") != -1)
    {
      concentration_units = line[headers.at("concentration_units")];
    }
    aqm.setConcentrationUnits(concentration_units);

    // statistics
    Int n_points = 0;
    if (headers.at("n_points") != -1)
    {
      n_points = (line[headers.at("n_points")].empty()) ? 0.0 : std::stoi(line[headers.at("n_points")]);
    }
    aqm.setNPoints(n_points);
    double correlation_coefficient = 0.0;
    if (headers.at("correlation_coefficient") != -1)
    {
      correlation_coefficient = (line[headers.at("correlation_coefficient")].empty()) ? 0.0 : std::stod(line[headers.at("correlation_coefficient")]);
    }
    aqm.setCorrelationCoefficient(correlation_coefficient);

    // transformation model
    String transformation_model = "";
    if (headers.at("transformation_model") != -1)
    {
      transformation_model = line[headers.at("transformation_model")];
    }
    Param transformation_model_params;
    for (auto const& kv : params_headers)
    {
      // cast doubles
      std::vector<String> param_doubles {"slope", "intercept", "wavelength", "span", "delta", "x_datum_min", "y_datum_min", "x_datum_max", "y_datum_max"}; 
      // cast integers
      std::vector<String> param_ints {"num_nodes", "boundary_condition", "num_iterations"};

      if (std::find(param_doubles.begin(), param_doubles.end(), kv.first) != param_doubles.end())
      {
        transformation_model_params.setValue(kv.first,std::stod(line[kv.second]));
      }      
      else if (std::find(param_ints.begin(), param_ints.end(), kv.first) != param_ints.end())
      {
        transformation_model_params.setValue(kv.first,std::stoi(line[kv.second]));
      }
      else
      {
        transformation_model_params.setValue(kv.first,line[kv.second]);
      }
      
    }
    aqm.setTransformationModel(transformation_model);
    aqm.setTransformationModelParams(transformation_model_params);
  }

  void AbsoluteQuantitationMethodFile::store(
    const String& filename,
    const std::vector<AbsoluteQuantitationMethod>& aqm_list
  )
  {
    clear(); // clear the buffer_
    const String headers = "IS_name,component_name,feature_name,concentration_units,llod,ulod,lloq,uloq,correlation_coefficient,n_points,transformation_model";
    StringList split_headers;
    headers.split(',', split_headers);
    StringList tm_params_names; // transformation model params
    if (aqm_list.size())
    {
      const Param tm_params = aqm_list[0].getTransformationModelParams();
      for (const Param::ParamEntry& param : tm_params)
      {
        if (std::find(tm_params_names.begin(), tm_params_names.end(), param.name) == tm_params_names.end())
        {
          tm_params_names.insert(tm_params_names.begin(), param.name);
        }
      }
      for (const String name : tm_params_names)
      {
        split_headers.push_back("transformation_model_param_" + name);
      }
    }
    addRow(split_headers);
    for (const AbsoluteQuantitationMethod& aqm : aqm_list)
    {
      StringList row(split_headers.size());
      row[0] = aqm.getISName();
      row[1] = aqm.getComponentName();
      row[2] = aqm.getFeatureName();
      row[3] = aqm.getConcentrationUnits();
      row[4] = aqm.getLLOD();
      row[5] = aqm.getULOD();
      row[6] = aqm.getLLOQ();
      row[7] = aqm.getULOQ();
      row[8] = aqm.getCorrelationCoefficient();
      row[9] = aqm.getNPoints();
      row[10] = aqm.getTransformationModel();
      const Param tm_params = aqm.getTransformationModelParams();
      for (Size i = 0, j = 11; i < tm_params_names.size(); ++i, ++j)
      {
        row[j] = tm_params.getValue(tm_params_names[i]);
      }
      addRow(row);
    }
    CsvFile::store(filename);
  }
} // namespace OpenMS
