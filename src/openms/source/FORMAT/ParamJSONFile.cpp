// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamJSONFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#if defined(ENABLE_TDL)
#include <tdl/tdl.h>
#else
#include <stdexcept>
#endif

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

// Reports a JSON value whose shape does not fit a parameter tagged as a file.
// This method is 'static' so no external linkage occurs
[[noreturn]] static void raiseFileShapeError(const std::string& key, const std::string& accepted, const json& node)
{
  std::string msg = "Parameter '" + key + "' is tagged as a file, but its JSON value is of type '" + std::string(node.type_name())
                    + "'. Accepted are " + accepted + ".";
  throw OpenMS::Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
}

// Extracts a single file (or directory) path from a JSON value.
// Accepted are a CWL 'File'/'Directory' object ({"class": "File", "path": "in.mzML"}), which is what a
// CWL runner writes, and a plain path string, which is what hand-written JSON and 'is_executable'
// parameters (they may name a tool on PATH rather than a file) use.
// This method is 'static' so no external linkage occurs
static std::string extractFilePath(const json& node, const std::string& key)
{
  if (node.is_string())
  {
    return node.get<std::string>();
  }
  if (node.is_object())
  {
    auto path = node.find("path");
    if (path != node.end() && path->is_string())
    {
      return path->get<std::string>();
    }
  }
  raiseFileShapeError(key, "a path string and an object with a 'path' string", node);
}

// Extracts a list of file (or directory) paths from a JSON value.
// A CWL runner serializes a 'File[]' input as an array of 'File' objects, so all of the following
// shapes have to be understood (see issue #10121):
//   (a) an array of CWL 'File'/'Directory' objects, written by every CWL runner:
//       [{"class": "File", "path": "a.mzML"}, {"class": "File", "path": "b.mzML"}]
//   (b) an array of plain path strings: ["a.mzML", "b.mzML"]
//   (c) a single object whose 'path' holds the array, which is what ParamJSONFile::store writes:
//       {"class": "File", "path": ["a.mzML", "b.mzML"]}
// This method is 'static' so no external linkage occurs
static std::vector<std::string> extractFilePathList(const json& node, const std::string& key)
{
  std::vector<std::string> paths;
  if (node.is_array()) // (a) and (b)
  {
    paths.reserve(node.size());
    for (const auto& element : node)
    {
      paths.push_back(extractFilePath(element, key));
    }
    return paths;
  }
  if (node.is_object()) // (c)
  {
    auto path = node.find("path");
    if (path != node.end())
    {
      if (path->is_array())
      {
        paths.reserve(path->size());
        for (const auto& element : *path)
        {
          paths.push_back(extractFilePath(element, key));
        }
        return paths;
      }
      if (path->is_string()) // a single file for a list-valued parameter
      {
        return {path->get<std::string>()};
      }
    }
  }
  raiseFileShapeError(key, "an array of path strings, an array of objects with a 'path' string, and an object with a 'path' array", node);
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
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else if (!File::readable(filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
    }
    try
    {
      json jsonNode = json::parse(std::ifstream {filename});
      auto traverseJSONTree = std::function<void(std::string currentKey, json& node)>{};
      traverseJSONTree = [&](std::string currentKey, json& node) {
        if (!node.is_object())
        {
          std::string msg = "Ignoring JSON file '" + filename + "' because of unexpected data type. Expecting a dictionary as type.";
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
        }

        for (const auto& child : node.items())
        {
          auto key = currentKey + replaceAll(child.key(), "__", ":"); // This converts __ to ':', but ':' would also be an accepted delimiter

          auto node = child.value();
          if (node.is_null())
          {
            continue; // No value given
          }
          // An object is a nested parameter section, unless it is the value of a parameter. Only
          // leaves exist as entries, so the parameter tree answers this for every shape a value
          // can take, the {"class": "File", "path": ...} objects a CWL runner passes included.
          if (node.is_object() && !param.exists(key)) {
            traverseJSONTree(key + ":", node);
            continue;
          }
          if (!param.exists(key))
          {
            std::string msg = "Parameter " + key + " passed to '" + traces.front().name + "' is invalid. To prevent usage of wrong defaults, please update/fix the parameters!";
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", msg);
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
            else if (entry.tags.contains("input file"))
            {
              // A CWL runner writes a 'class: File' object, hand-written JSON a plain string. An
              // 'is_executable' parameter uses a string for a tool that is looked up on PATH.
              value = extractFilePath(node, key);
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
            if (entry.tags.contains("input file") || entry.tags.contains("output file"))
            {
              value = extractFilePathList(node, key);
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
      };
      traverseJSONTree(toolNamespace, jsonNode);
    }
    catch (const json::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", e.what());
      return false;
    }
    return true;
  }

  void ParamJSONFile::store(const std::string& filename, const Param& param, const ToolInfo&) const
  {
    std::ofstream os;
    std::ostream* os_ptr;
    if (filename != "-")
    {
      os.open(filename.c_str(), std::ofstream::out);
      if (!os)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      os_ptr = &os;
    }
    else
    {
      os_ptr = &std::cout;
    }

    writeToStream(os_ptr, param);
  }

  void ParamJSONFile::writeToStream(std::ostream* os_ptr, const Param& param) const
  {
#if defined(ENABLE_TDL)
    std::ostream& os = *os_ptr;

    // discover the name of the first nesting Level
    // this is expected to result in something like "ToolName:1:"
    auto traces = param.begin().getTrace();
    std::string toolNamespace = traces.front().name + ":1:";

    std::vector<tdl::Node> stack;
    stack.push_back(tdl::Node {});

    json jsonDoc{};


    auto param_it = param.begin();
    for (auto last = param.end(); param_it != last; ++param_it)
    {
      for (auto& trace : param_it.getTrace())
      {
        if (trace.opened)
        {
          stack.push_back(tdl::Node {trace.name, trace.description, {}, tdl::Node::Children {}});
        }
        else // these nodes must be closed
        {
          auto top = stack.back();
          stack.pop_back();
          auto& children = std::get<tdl::Node::Children>(stack.back().value);
          children.push_back(top);
        }
      }

      // converting trags to tdl compatible tags
      std::set<std::string> tags;
      for (auto const& t : param_it->tags)
      {
        if (t == "input file")
        {
          tags.insert("file");
        }
        else if (t == "output file")
        {
          tags.insert("file");
          tags.insert("output");
        }
        else if (t == "output prefix")
        {
          tags.insert("output");
          tags.insert("prefixed");
        }
        else
        {
          tags.insert(t);
        }
      }

      // Sets a single value into the tdl library
      auto name = param_it->name;
      if (stack.size() > 2) {
        json node{};

        switch (param_it->value.valueType())
        {
          case ParamValue::INT_VALUE:
            node = static_cast<int64_t>(param_it->value);
            break;
          case ParamValue::DOUBLE_VALUE:
            node = static_cast<double>(param_it->value);
            break;
          case ParamValue::STRING_VALUE:
            if ((param_it->valid_strings.size() == 2 && param_it->valid_strings[0] == "true" && param_it->valid_strings[1] == "false")
               || (param_it->valid_strings.size() == 2 && param_it->valid_strings[0] == "false" && param_it->valid_strings[1] == "true"))
            {
                node = param_it->value.toBool();
            } else {
                if (tags.contains("file") && !tags.contains("output")) {
                    node["class"] = "File";
                    node["path"] = param_it->value.toString();
                } else {
                    node = param_it->value.toString();
                }
            }
            break;
          case ParamValue::INT_LIST:
            node = param_it->value.toIntVector();
            break;
          case ParamValue::DOUBLE_LIST:
            node = param_it->value.toDoubleVector();
            break;
          case ParamValue::STRING_LIST:
            if (tags.contains("file") && !tags.contains("output")) {
                node["class"] = "File";
                node["path"] = param_it->value.toStringVector();
            } else {
                node = param_it->value.toStringVector();
            }
            break;
          default:
            break;
        }

        // Add newly created node to json document
        if (!flatHierarchy) {
          // Traverse to the correct node
          auto* parent = &jsonDoc;
          for (size_t i{3}; i < stack.size(); ++i) {
              parent = &(*parent)[stack[i].name];
          }
          (*parent)[name] = node;
        } else {
          // Expand name to include all namespaces
          for (size_t i{0}; i < stack.size()-3; ++i) {
            auto const& e = stack[stack.size()-1-i];
            name = e.name + "__" + name;
          }
          jsonDoc[name] = node;
        }
      }
    }

    while (stack.size() > 1)
    {
      auto top = stack.back();
      stack.pop_back();
      auto& children = std::get<tdl::Node::Children>(stack.back().value);
      children.push_back(top);
    }
    assert(stack.size() == 1);

    os << jsonDoc.dump(2);
#else
    (void)os_ptr;
    (void)param;
    throw std::runtime_error{"TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature."};
#endif
  }
} // namespace OpenMS
