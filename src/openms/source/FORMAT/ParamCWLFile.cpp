// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamCWLFile.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>

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
  void ParamCWLFile::load(const std::string& filename, Param& param)
  {
    // discover the name of the first nesting Level
    // this is expected to result in something like "ToolName:1:"
    auto traces = param.begin().getTrace();
    std::string toolNamespace = traces.front().name + ":1:";

    json jsonNode = json::parse(std::ifstream {filename});
    if (!jsonNode.is_object())
    {
      std::string msg = "Ignoring JSON file '" + filename + "' because of unexpected data type. Expecting a dictionary as root type.";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
    }
    for (auto child : jsonNode.items())
    {
      auto key = child.key();
      key = toolNamespace + replaceAll(key, "__", ":"); // This converts __ to ':', but ':' would also be an accepted delimiter

      auto node = child.value();
      if (node.is_null())
        continue; // No value given
      if (!param.exists(key))
      {
        OPENMS_LOG_WARN << "Unknown (or deprecated) Parameter '" << key << "' given in outdated parameter file! Ignoring parameter." << std::endl;
        continue;
      }
      auto const& entry = param.getEntry(key);
      auto value = entry.value;
      if (entry.value.valueType() == ParamValue::ValueType::STRING_VALUE)
      {
        if ((entry.valid_strings.size() == 2 && entry.valid_strings[0] == "true" && entry.valid_strings[1] == "false")
           || (entry.valid_strings.size() == 2 && entry.valid_strings[0] == "false" && entry.valid_strings[1] == "true"))
        {
          value = node.get<bool>() ? "true" : "false";
        }
        else if (entry.tags.count("input file") || entry.tags.count("output file"))
        {
          value = node["path"].get<std::string>();
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
        if (entry.tags.count("input file") || entry.tags.count("output file"))
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

} // namespace OpenMS
