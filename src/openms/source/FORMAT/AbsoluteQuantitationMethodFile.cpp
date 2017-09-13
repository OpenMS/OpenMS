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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

namespace OpenMS
{

  AbsoluteQuantitationMethodFile::AbsoluteQuantitationMethodFile()
  {
  }

  AbsoluteQuantitationMethodFile::~AbsoluteQuantitationMethodFile()
  {
  }

  void AbsoluteQuantitationMethodFile::load(const String & filename, AbsoluteQuantitationMethod & aqm)
  {
    // read in the .csv file
    char is = ',';
    bool ie = false; 
    Int first_n = -1;
    fload(filename, is, ie, first_n);

    // parse the file
    std::map<std::string,int> headers;
    std::map<std::string,int> params_headers;
    StringList line, header;
    for (size_t i = 0; i < CsvFile::rowCount(); ++i)
    {
      if (i == 0) // header row
      {
        CsvFile::getRow(i, header);
        parseHeader(header, headers, params_headers);
      }
      else
      {
        CsvFile::getRow(i, line);
        parseLine(line, headers, params_headers, aqm);    
      }      
    }
  }

  void AbsoluteQuantitationMethodFile::parseHeader(StringList & line, std::map<std::string, int> & headers,
    std::map<std::string, int> & params_headers)
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
    headers["actual_concentration"] = -1;
    headers["n_points"] = -1;
    headers["transformation_model"] = -1;
    std::string param_header = "transformation_model_param_";
    
    // parse the header columns
    for (size_t i = 0; i < line.size(); ++i)
    {
      // parse transformation_model_params
      if (line[i].find(param_header) != std::string::npos) 
      {
        line[i].erase(line[i].begin()+param_header.size()); 
        params_headers[line[i]] = i;
        std::cout << line[i] << std::endl;
      }      
      else // parse all other header entries
      {
        headers[line[i]] = i;
      }
    }

    // // default headers
    // std::map<std::string,std::string> string_map;
    // string_map["IS_name"] = "";
    // string_map["component_name"] = "";
    // string_map["feature_name"] = "";
    // string_map["concentration_units"] = "";
    // string_map["transformation_model"] = "";
    // std::map<std::string,double> float_map;
    // float_map["llod"] = 0.0;
    // float_map["ulod"] = 0.0;
    // float_map["lloq"] = 0.0;
    // float_map["uloq"] = 0.0;
    // float_map["correlation_coefficient"] = 0.0;
    // float_map["actual_concentration"] = 0.0;
    // std::map<std::string,int> int_map;
    // int_map["n_points"] = 0;

  }

  void AbsoluteQuantitationMethodFile::parseLine(StringList & line, std::map<std::string,int> & headers, 
    std::map<std::string,int> & params_headers, AbsoluteQuantitationMethod & aqm)
  {
    // component, IS, and feature names
    std::string component_name = "";
    if (headers["component_name"] != -1)
    {
      component_name = line[headers["component_name"]];
    }
    std::string feature_name = "";
    if (headers["feature_name"] != -1)
    {
      feature_name = line[headers["feature_name"]];
    }
    std::string IS_name = "";
    if (headers["IS_name"] != -1)
    {
      IS_name = line[headers["IS_name"]];
    }
    aqm.setComponentISFeatureNames(component_name, IS_name, feature_name);

    // LODs
    double llod = 0.0;
    if (headers["llod"] != -1)
    {
      llod = std::stod(line[headers["llod"]]);
    }
    double ulod = 0.0;
    if (headers["ulod"] != -1)
    {
      ulod = std::stod(line[headers["ulod"]]);
    }
    aqm.setLOD(llod,ulod);

    // LOQs
    double lloq = 0.0;
    if (headers["lloq"] != -1)
    {
      lloq = std::stod(line[headers["lloq"]]);
    }
    double uloq = 0.0;
    if (headers["uloq"] != -1)
    {
      uloq = std::stod(line[headers["uloq"]]);
    }
    aqm.setLOQ(lloq,uloq);

    // actual concentration
    double actual_concentration = 0.0;
    if (headers["actual_concentration"] != -1)
    {
      actual_concentration = std::stod(line[headers["actual_concentration"]]);
    }
    aqm.setActualConcentration(actual_concentration);

    // concentration units
    std::string concentration_units = "";
    if (headers["concentration_units"] != -1)
    {
      concentration_units = line[headers["concentration_units"]];
    }
    aqm.setConcentrationUnits(concentration_units);

    // statistics
    int n_points = 0;
    if (headers["n_points"] != -1)
    {
      n_points = std::stoi(line[headers["n_points"]]);
    }
    double correlation_coefficient = 0.0;
    if (headers["correlation_coefficient"] != -1)
    {
      correlation_coefficient = std::stod(line[headers["correlation_coefficient"]]);
    }
    aqm.setStatistics(n_points, correlation_coefficient);

    // transformation model
    std::string transformation_model = "";
    if (headers["transformation_model"] != -1)
    {
      transformation_model = line[headers["transformation_model"]];
    }
    Param transformation_model_params;
    for (auto const& kv : params_headers)
    {
      // cast doubles
      std::vector<std::string> param_doubles {"slope", "intercept", "wavelength", "span", "delta", "x_datum_min", "y_datum_min", "x_datum_max", "y_datum_max"};      
      if (std::find(param_doubles.begin(), param_doubles.end(), kv.first) != param_doubles.end())
      {
        transformation_model_params.setValue(kv.first,std::stod(line[kv.second]));
      }
      // cast integers
      std::vector<std::string> param_ints {"num_nodes", "boundary_condition", "num_iterations"};      
      else if (std::find(param_ints.begin(), param_ints.end(), kv.first) != param_ints.end())
      {
        transformation_model_params.setValue(kv.first,std::stoi(line[kv.second]));
      }
      else
      {
        transformation_model_params.setValue(kv.first,line[kv.second]);
      }
      
    }
    aqm.setTransformationModel(transformation_model, transformation_model_params);
  }

  void AbsoluteQuantitationMethodFile::store(const String & filename, const AbsoluteQuantitationMethod & aqm)
  {
    // TODO
  }

} // namespace OpenMS
