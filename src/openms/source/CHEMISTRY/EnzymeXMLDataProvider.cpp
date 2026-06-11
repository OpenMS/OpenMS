// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymeXMLDataProvider.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>

namespace OpenMS
{
  template <typename EnzymeType>
  EnzymeXMLDataProvider<EnzymeType>::EnzymeXMLDataProvider(const std::string& filename, bool optional)
    : filename_(filename), optional_(optional)
  {
  }

  template <typename EnzymeType>
  std::vector<std::unique_ptr<EnzymeType>> EnzymeXMLDataProvider<EnzymeType>::loadEnzymes()
  {
    std::vector<std::unique_ptr<EnzymeType>> enzymes;

    std::string file;
    try
    {
      file = File::find(filename_);
    }
    catch (Exception::FileNotFound&)
    {
      if (optional_)
      {
        return enzymes; // empty
      }
      throw;
    }

    Param param;
    ParamXMLFile().load(file, param);
    if (param.empty())
    {
      return enzymes;
    }

    std::vector<std::string> split;
    StringUtils::split(std::string(param.begin().getName()), ':', split);
    if (split[0] != "Enzymes")
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, split[0], "name 'Enzymes' expected");
    }

    try
    {
      // helper: construct an enzyme from accumulated key/value pairs
      auto build_enzyme = [&enzymes](const std::map<std::string, std::string>& kv)
      {
        auto enzyme = std::make_unique<EnzymeType>();
        for (const auto& [key, value] : kv)
        {
          if (!enzyme->setValueFromFile(key, value))
          {
            OPENMS_LOG_ERROR << "Error while parsing enzymes file: unknown key '" << key << "' with value '" << value << "'" << std::endl;
          }
        }
        enzymes.push_back(std::move(enzyme));
      };

      std::map<std::string, std::string> values;
      std::string previous_enzyme = split[1];
      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        StringUtils::split(std::string(it.getName()), ':', split);
        if (split[0] != "Enzymes")
        {
          break;
        }
        if (split[1] != previous_enzyme)
        {
          build_enzyme(values);
          previous_enzyme = split[1];
          values.clear();
        }
        values[it.getName()] =std::string(it->value.toString());
      }
      build_enzyme(values); // add last enzyme
    }
    catch (Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), "");
    }

    return enzymes;
  }

  // Explicit template instantiations
  template class EnzymeXMLDataProvider<DigestionEnzymeProtein>;
  template class EnzymeXMLDataProvider<DigestionEnzymeRNA>;

} // namespace OpenMS
