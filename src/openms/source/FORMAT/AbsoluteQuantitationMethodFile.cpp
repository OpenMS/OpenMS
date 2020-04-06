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

#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>
#include <fstream>
#include <boost/regex.hpp>

namespace OpenMS
{
  void AbsoluteQuantitationMethodFile::load(const String & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list)
  {
    aqm_list.clear();
    CsvFile::load(filename, ',', false, -1);
    std::map<String, Size> headers;
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
    const std::map<String, Size> & headers,
    AbsoluteQuantitationMethod & aqm
  ) const
  {
    StringList tl = line; // trimmed line
    for (String& s : tl)
    {
      s.trim();
    }
    aqm.setComponentName(headers.count("component_name") ? tl[headers.at("component_name")] : "");
    aqm.setFeatureName(headers.count("feature_name") ? tl[headers.at("feature_name")] : "");
    aqm.setISName(headers.count("IS_name") ? tl[headers.at("IS_name")] : "");
    aqm.setLLOD(!headers.count("llod") || tl[headers.at("llod")].empty() ? 0 : std::stod(tl[headers.at("llod")]));
    aqm.setULOD(!headers.count("ulod") || tl[headers.at("ulod")].empty() ? 0 : std::stod(tl[headers.at("ulod")]));
    aqm.setLLOQ(!headers.count("lloq") || tl[headers.at("lloq")].empty() ? 0 : std::stod(tl[headers.at("lloq")]));
    aqm.setULOQ(!headers.count("uloq") || tl[headers.at("uloq")].empty() ? 0 : std::stod(tl[headers.at("uloq")]));
    aqm.setConcentrationUnits(headers.count("concentration_units") ? tl[headers.at("concentration_units")] : "");
    aqm.setNPoints(!headers.count("n_points") || tl[headers.at("n_points")].empty() ? 0 : std::stoi(tl[headers.at("n_points")]));
    aqm.setCorrelationCoefficient(
      !headers.count("correlation_coefficient") || tl[headers.at("correlation_coefficient")].empty()
        ? 0
        : std::stod(tl[headers.at("correlation_coefficient")])
    );
    aqm.setTransformationModel(headers.count("transformation_model") ? tl[headers.at("transformation_model")] : "");
    Param tm_params;
    for (const std::pair<const String, Size>& h : headers)
    {
      const String& header = h.first;
      const Size& i = h.second;
      boost::smatch m;
      if (boost::regex_search(header, m, boost::regex("transformation_model_param_(.+)")))
      {
        setCastValue_(String(m[1]), tl[i], tm_params);
      }
    }
    aqm.setTransformationModelParams(tm_params);
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
        row[j] = tm_params.exists(tm_params_names[i]) ? tm_params.getValue(tm_params_names[i]) : "";
      }
      addRow(row);
    }
    CsvFile::store(filename);
  }

  void AbsoluteQuantitationMethodFile::setCastValue_(const String& key, const String& value, Param& params) const
  {
    const std::vector<String> param_doubles {
      "slope", "intercept", "wavelength", "span", "delta", "x_datum_min", "y_datum_min", "x_datum_max", "y_datum_max"
    };
    const std::vector<String> param_ints {"num_nodes", "boundary_condition", "num_iterations"};
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
