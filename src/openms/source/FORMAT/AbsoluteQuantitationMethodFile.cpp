// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>
#include <boost/regex.hpp>

namespace OpenMS
{
  void AbsoluteQuantitationMethodFile::load(const std::string & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)
  {
    aqm_list.clear();
    CsvFile::load(filename, ',', false, -1);
    std::map<std::string, Size> headers;
    StringList sl;
    if (rowCount() >= 2) // no need to read headers if that's the only line inside the file
    {
      getRow(0, sl);
      for (Size i = 0; i < sl.size(); ++i)
      {
        headers[sl[i]] = i; // for each header found, assign an index value to it
      }
      if (!( // if any of these headers is missing, warn the user
        headers.count("IS_name") &&
        headers.count("component_name") &&
        headers.count("feature_name") &&
        headers.count("concentration_units") &&
        headers.count("llod") &&
        headers.count("ulod") &&
        headers.count("lloq") &&
        headers.count("uloq") &&
        headers.count("correlation_coefficient") &&
        headers.count("n_points") &&
        headers.count("transformation_model")
      ))
      {
        OPENMS_LOG_WARN << "One or more of the following columns are missing:\n";
        OPENMS_LOG_WARN << "IS_name\n";
        OPENMS_LOG_WARN << "component_name\n";
        OPENMS_LOG_WARN << "feature_name\n";
        OPENMS_LOG_WARN << "concentration_units\n";
        OPENMS_LOG_WARN << "llod\n";
        OPENMS_LOG_WARN << "ulod\n";
        OPENMS_LOG_WARN << "lloq\n";
        OPENMS_LOG_WARN << "uloq\n";
        OPENMS_LOG_WARN << "correlation_coefficient\n";
        OPENMS_LOG_WARN << "n_points\n";
        OPENMS_LOG_WARN << "transformation_model\n" << std::endl;
      }
    }
    for (Size i = 1; i < rowCount(); ++i)
    {
      getRow(i, sl);
      AbsoluteQuantitationMethod aqm;
      parseLine_(sl, headers, aqm);
      aqm_list.push_back(aqm);
    }
  }

  void AbsoluteQuantitationMethodFile::parseLine_(
    const StringList & line,
    const std::map<std::string, Size> & headers,
    AbsoluteQuantitationMethod & aqm
  ) const
  {
    StringList tl = line; // trimmed line
    for (std::string& s : tl)
    {
      StringUtils::trim(s);
    }
    aqm.setComponentName(headers.count("component_name") ? tl[headers.at("component_name")] : "");
    aqm.setFeatureName(headers.count("feature_name") ? tl[headers.at("feature_name")] : "");
    aqm.setISName(headers.count("IS_name") ? tl[headers.at("IS_name")] : "");
    aqm.setLLOD(!headers.contains("llod") || tl[headers.at("llod")].empty() ? 0 : std::stod(tl[headers.at("llod")]));
    aqm.setULOD(!headers.contains("ulod") || tl[headers.at("ulod")].empty() ? 0 : std::stod(tl[headers.at("ulod")]));
    aqm.setLLOQ(!headers.contains("lloq") || tl[headers.at("lloq")].empty() ? 0 : std::stod(tl[headers.at("lloq")]));
    aqm.setULOQ(!headers.contains("uloq") || tl[headers.at("uloq")].empty() ? 0 : std::stod(tl[headers.at("uloq")]));
    aqm.setConcentrationUnits(headers.count("concentration_units") ? tl[headers.at("concentration_units")] : "");
    aqm.setNPoints(!headers.contains("n_points") || tl[headers.at("n_points")].empty() ? 0 : std::stoi(tl[headers.at("n_points")]));
    aqm.setCorrelationCoefficient(
      !headers.contains("correlation_coefficient") || tl[headers.at("correlation_coefficient")].empty()
        ? 0
        : std::stod(tl[headers.at("correlation_coefficient")])
    );
    aqm.setTransformationModel(headers.count("transformation_model") ? tl[headers.at("transformation_model")] : "");
    Param tm_params;
    for (const std::pair<const std::string, Size>& h : headers)
    {
      const std::string& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("transformation_model_param_(.+)")))
      {
        setCastValue_(StringUtils::toStr(m[1]), tl[i], tm_params);
      }
    }
    aqm.setTransformationModelParams(tm_params);
  }

  void AbsoluteQuantitationMethodFile::store(
    const std::string& filename,
    const std::vector<AbsoluteQuantitationMethod>& aqm_list
  )
  {
    clear(); // clear the buffer_
    const std::string headers = "IS_name,component_name,feature_name,concentration_units,llod,ulod,lloq,uloq,correlation_coefficient,n_points,transformation_model";
    StringList split_headers;
    StringUtils::split(headers, ',', split_headers);
    StringList tm_params_names; // transformation model params
    if (!aqm_list.empty())
    {
      const Param tm_params = aqm_list[0].getTransformationModelParams();
      for (const Param::ParamEntry& param : tm_params)
      {
        tm_params_names.insert(tm_params_names.begin(), param.name);
        split_headers.insert(split_headers.begin() + 11, "transformation_model_param_" + param.name);
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
        row[j] = tm_params.exists(tm_params_names[i]) ? tm_params.getValue(tm_params_names[i]).toString() : "";
      }
      addRow(row);
    }
    CsvFile::store(filename);
  }

  void AbsoluteQuantitationMethodFile::setCastValue_(const std::string& key, const std::string& value, Param& params) const
  {
    const std::vector<std::string> param_doubles {
      "slope", "intercept", "wavelength", "span", "delta", "x_datum_min", "y_datum_min", "x_datum_max", "y_datum_max"
    };
    const std::vector<std::string> param_ints {"num_nodes", "boundary_condition", "num_iterations"};
    if (std::find(param_doubles.begin(), param_doubles.end(), key) != param_doubles.end())
    {
      params.setValue(key, value.empty() ? 0 : std::stod(value));
    }
    else if (std::find(param_ints.begin(), param_ints.end(), key) != param_ints.end())
    {
      params.setValue(key, value.empty() ? 0 : std::stoi(value));
    }
    else
    {
      params.setValue(key,value);
    }
  }
} // namespace OpenMS
