// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamCWLFile.h>
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
  void ParamCWLFile::store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const
  {
    std::ofstream os;
    std::ostream* os_ptr;
    if (filename != "-")
    {
      os.open(filename.c_str(), std::ofstream::out);
      if (!os)
      {
        // Replace the OpenMS specific exception with a std exception
        // Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
        throw std::ios::failure("Unable to create file: " + filename);
      }
      os_ptr = &os;
    }
    else
    {
      os_ptr = &std::cout;
    }

    writeCWLToStream(os_ptr, param, tool_info);
  }

  void ParamCWLFile::writeCWLToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const
  {
    std::ostream& os = *os_ptr;
    os.precision(std::numeric_limits<double>::digits10);

    tdl::ToolInfo tdl_tool_info;
    tdl_tool_info.metaInfo.version = tool_info.version_;
    tdl_tool_info.metaInfo.name = tool_info.name_;
    tdl_tool_info.metaInfo.docurl = tool_info.docurl_;
    tdl_tool_info.metaInfo.category = tool_info.category_;
    tdl_tool_info.metaInfo.description = tool_info.description_;
    for (auto cite : tool_info.citations_)
    {
      tdl::Citation tdl_citation;
      tdl_citation.doi = cite;
      tdl_citation.url = "";
      tdl_tool_info.metaInfo.citations.push_back(tdl_citation);
    }

    // discover the name of the first nesting Level
    // this is expected to result in something like "ToolName:1:"
    auto traces = param.begin().getTrace();
    std::string toolNamespace = traces.front().name + ":1:";

    std::vector<tdl::Node> stack;
    stack.push_back(tdl::Node {});

    auto param_it = param.begin();
    for (auto last = param.end(); param_it != last; ++param_it)
    {
      for (auto& trace : param_it.getTrace())
      {
        if (trace.opened)
        {
          // First nested param should be the executable name of the tool
          if (tdl_tool_info.metaInfo.executableName.empty())
          {
            tdl_tool_info.metaInfo.executableName = trace.name;
          }
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

      // converting OpenMS tags to tdl compatible tags
      std::set<std::string> tags;
      for (auto const& t : param_it->tags)
      {
        if (t == TOPPBase::TAG_INPUT_FILE)
        {
          tags.insert("file");
        }
        else if (t == TOPPBase::TAG_OUTPUT_FILE)
        {
          tags.insert("file");
          tags.insert("output");
        }
        else if (t == TOPPBase::TAG_OUTPUT_PREFIX)
        {
          tags.insert("output");
          tags.insert("prefixed");
        }
        else if (t == TOPPBase::TAG_OUTPUT_DIR)
        {
          tags.insert("directory");
          tags.insert("output");
        }
        else
        {
          tags.insert(t);
        }
      }

      // Sets a single value into the tdl library
      auto setValue = [&](auto value) { std::get<tdl::Node::Children>(stack.back().value).push_back(tdl::Node {param_it->name, param_it->description, tags, value}); };
      // Sets a value including their limits into the tdl library
      auto setValueLimits = [&](auto value) {
        using T = typename decltype(value.minLimit)::value_type;
        if (value.minLimit == -std::numeric_limits<T>::max())
          value.minLimit.reset();
        if (value.maxLimit == std::numeric_limits<T>::max())
          value.maxLimit.reset();
        setValue(value);
      };
      switch (param_it->value.valueType())
      {
        case ParamValue::INT_VALUE:
          setValueLimits(tdl::IntValue {static_cast<int>(param_it->value), param_it->min_int, param_it->max_int});
          break;
        case ParamValue::DOUBLE_VALUE:
          setValueLimits(tdl::DoubleValue {static_cast<double>(param_it->value), param_it->min_float, param_it->max_float});
          break;
        case ParamValue::STRING_VALUE:
          if (param_it->valid_strings.size() == 2 && param_it->valid_strings[0] == "true" && param_it->valid_strings[1] == "false" && param_it->value == "false")
          {
            std::get<tdl::Node::Children>(stack.back().value).push_back(tdl::Node {param_it->name, param_it->description, tags, false});
          }
          else
          {
            std::get<tdl::Node::Children>(stack.back().value)
              .push_back(tdl::Node {param_it->name, param_it->description, tags, tdl::StringValue {static_cast<std::string>(param_it->value), param_it->valid_strings}});
          }
          break;
        case ParamValue::INT_LIST: {
          auto value = tdl::IntValueList {};
          value.value = param_it->value.toIntVector();
          value.minLimit = param_it->min_int;
          value.maxLimit = param_it->max_int;
          setValueLimits(value);
          break;
        }
        case ParamValue::DOUBLE_LIST: {
          auto value = tdl::DoubleValueList {};
          value.value = param_it->value.toDoubleVector();
          value.minLimit = param_it->min_float;
          value.maxLimit = param_it->max_float;
          setValueLimits(value);
          break;
        }
        case ParamValue::STRING_LIST: {
          auto value = tdl::StringValueList {};
          value.value = param_it->value.toStringVector();
          value.validValues = param_it->valid_strings;
          auto param = tdl::Node {param_it->name, param_it->description, tags, value};
          std::get<tdl::Node::Children>(stack.back().value).push_back(param);
          break;
        }
        default:
          break;
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

    auto allChildren = tdl::Node::Children{};
    // fix naming of all children, by appending their parents name
    // skipping the first two levels, since they are always the "ToolName:1:" keys
    auto renameNodes = std::function<void(tdl::Node&, tdl::Node const*, int)> {};
    renameNodes = [&](tdl::Node& element, tdl::Node const* parent, int level) {
      if (parent != nullptr && !parent->name.empty())
      {
        element.name = parent->name + ":" + element.name;
      }

      if (auto children = std::get_if<tdl::Node::Children>(&element.value))
      {
        for (auto& child : *children)
        {
          renameNodes(child, &element, level + 1);
        }
      } else if (!element.name.empty()) {
        allChildren.push_back(element);
      }
    };
    renameNodes(stack.back(), nullptr, 0);

    if (flatHierarchy) {
        stack = allChildren;
    }

    // This does different things
    // 1. uses a safer sign than ':' for output
    // 2. strips of the toolNamespace
    // 3. ignore certain options
    //
    // returns true if object can be removed
    auto recursiveCleanup = std::function<void(tdl::Node&)> {};
    recursiveCleanup = [&](tdl::Node& element) {
      auto& name = element.name;

      // strip of the tool namespace part and ignore entries that aren't part of the name space (like ToolName:version)
      if (name.size() >= toolNamespace.size() && name.substr(0, toolNamespace.size()) == toolNamespace)
      {
        name = name.substr(toolNamespace.size());
      }
      else
      {
        name = "";
      }

      if (flatHierarchy) {
        // replace all ':' with '__'
        name = replaceAll(name, ":", "__");
      } else {
        if (auto pos = name.rfind(':'); pos != std::string::npos) {
            name = name.substr(pos+1);
        }
      }

      // cleanup recursivly
      if (auto children = std::get_if<tdl::Node::Children>(&element.value))
      {
        for (auto& child : *children) {
            recursiveCleanup(child);
        }
      }
    };
    for (auto& s : stack) {
        recursiveCleanup(s);
    }

    // recursiveMarkingRemove
    auto recursiveMarkRemove = std::function<bool(tdl::Node&)> {};
    recursiveMarkRemove = [&](tdl::Node& element) -> bool {
      // remove children that are marked for cleanup
      if (auto children = std::get_if<tdl::Node::Children>(&element.value))
      {
        children->erase(std::remove_if(children->begin(), children->end(), [&](auto& child) {
                                            return recursiveMarkRemove(child);
                                        }),
                                        children->end());
         return children->empty() && element.name.empty();
      }
      return element.name.empty();
    };

    // Erase invalid entries
    stack.erase(std::remove_if(stack.begin(), stack.end(), [&](auto& child) {
                                        return recursiveMarkRemove(child);
                                    }),
                                    stack.end());


    // unroll nested options, until you have some with names
    while (stack.size() == 1 && stack.back().name.empty() && std::holds_alternative<tdl::Node::Children>(stack.back().value)) {
        auto v = std::get<tdl::Node::Children>(stack.back().value);
        stack = v;
    }

    tdl_tool_info.params = stack;

    // Removing the fake cli methods
    tdl::post_process_cwl = [&](YAML::Node node) {
      node["requirements"] = YAML::Load(R"-(
            InlineJavascriptRequirement: {}
            InitialWorkDirRequirement:
              listing:
                - entryname: cwl_inputs.json
                  entry: $(JSON.stringify(inputs))
        )-");
      node["arguments"] = YAML::Load(R"-(
            - -ini
            - cwl_inputs.json
        )-");
    };
    os << "# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin\n"
          "# SPDX-License-Identifier: Apache-2.0\n";

    os << convertToCWL(tdl_tool_info) << "\n";
  }
} // namespace OpenMS
