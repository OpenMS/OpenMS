// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamJSONFile.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <tdl/tdl.h>

using json = nlohmann::json;

// replaces the substrings inside a string
// This method is 'static' so no external linkage occurs
static std::string replaceAll(std::string str, const std::string& pattern, const std::string& replacement)
{
  size_t pos {0};
  while ((pos = str.find(pattern, pos)) != std::string::npos)
  {
    str.replace(pos, pattern.size(), replacement);
    pos += replacement.size();
  }
  return str;
}

namespace OpenMS
{
  bool ParamJSONFile::load(const std::string& filename, Param& param)
  {
    // discover the name of the first nesting Level
    // this is expected to result in something like "ToolName:1:"
    auto traces = param.begin().getTrace();
    std::string toolNamespace = traces.front().name + ":1:";

    std::ifstream ifs {filename};
    if (!ifs.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    try
    {
      json jsonNode = json::parse(std::ifstream {filename});
      if (!jsonNode.is_object())
      {
        std::string msg = "Ignoring JSON file '" + filename + "' because of unexpected data type. Expecting a dictionary as root type.";
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
      }
      for (const auto& child : jsonNode.items())
      {
        auto key = child.key();
        key = toolNamespace + replaceAll(key, "__", ":"); // This converts __ to ':', but ':' would also be an accepted delimiter

        auto node = child.value();
        if (node.is_null())
        {
          continue; // No value given
        }
        if (!param.exists(key))
        {
          OPENMS_LOG_ERROR << "Parameter " << key << " passed to '" << traces.front().name << "' is invalid. To prevent usage of wrong defaults, please update/fix the parameters!" << std::endl;
          return false;
        }
        auto const& entry = param.getEntry(key);
        auto value = entry.value;
        if (entry.value.valueType() == ParamValue::ValueType::STRING_VALUE)
        {
          if ((entry.valid_strings.size() == 2 && entry.valid_strings[0] == "true" && entry.valid_strings[1] == "false") ||
              (entry.valid_strings.size() == 2 && entry.valid_strings[0] == "false" && entry.valid_strings[1] == "true"))
          {
            value = node.get<bool>() ? "true" : "false";
          }
          else if (entry.tags.count("input file"))
          {
            // If this is an input file and 'is_executable' is set. this can be of 'class: File' or 'type: string'
            if (entry.tags.count("is_executable"))
            {
              if (node.is_object())
              {
                value = node["path"].get<std::string>();
              }
              else
              {
                value = node.get<std::string>();
              }
            }
            // Just a normal input file
            else
            {
              value = node["path"].get<std::string>();
            }
          }
          else
          {
            value = node.get<std::string>();
          }
        }
        else if (entry.value.valueType() == ParamValue::ValueType::INT_VALUE)
        {
          value = node.get<int64_t>();
        }
        else if (entry.value.valueType() == ParamValue::ValueType::DOUBLE_VALUE)
        {
          value = node.get<double>();
        }
        else if (entry.value.valueType() == ParamValue::ValueType::STRING_LIST)
        {
          if (entry.tags.count("input file"))
          {
            value = node["path"].get<std::vector<std::string>>();
          }
          else
          {
            value = node.get<std::vector<std::string>>();
          }
        }
        else if (entry.value.valueType() == ParamValue::ValueType::INT_LIST)
        {
          value = node.get<std::vector<int>>();
        }
        else if (entry.value.valueType() == ParamValue::ValueType::DOUBLE_LIST)
        {
          value = node.get<std::vector<double>>();
        }
        else if (entry.value.valueType() == ParamValue::ValueType::EMPTY_VALUE)
        {
          // Nothing happens here
          OPENMS_LOG_WARN << "Ignoring entry '" << key << "' because of unknown type 'EMPTY_VALUE'." << std::endl;
        }
        param.setValue(key, value);
      }
    }
    catch (const json::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", e.what());
    }
    return true;
  }
} // namespace OpenMS
