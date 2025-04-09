// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ruben Grünberg $
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParamCTDFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>


namespace OpenMS
{
  void ParamCTDFile::store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const
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

    writeCTDToStream(os_ptr, param, tool_info);
  }

  void ParamCTDFile::writeCTDToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const
  {
    std::ostream& os = *os_ptr;
    os.precision(std::numeric_limits<double>::digits10);

    // write ctd specific stuff
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    os << R"(<tool ctdVersion="1.8" version=")" << tool_info.version_ << R"(" name=")" << tool_info.name_ << R"(" docurl=")" << tool_info.docurl_ << R"(" category=")" << tool_info.category_
       << "\" >\n";
    os << "<description><![CDATA[" << tool_info.description_ << "]]></description>\n";
    os << "<manual><![CDATA[" << tool_info.description_ << "]]></manual>\n";
    os << "<citations>\n";

    for (auto& doi : tool_info.citations_)
    {
      os << "  <citation doi=\"" << doi << "\" url=\"\" />\n";
    }

    os << "</citations>\n";
    os << "<PARAMETERS version=\"" << schema_version_ << R"(" xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS)" << schema_location_
       << "\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // Write the xml stuff
    uint32_t indentations = 2;
    auto param_it = param.begin();
    for (auto last = param.end(); param_it != last; ++param_it)
    {
      for (auto& trace : param_it.getTrace())
      {
        if (trace.opened)
        {
          std::string d = trace.description;
          replace(d, '\n', "#br#");

          os << std::string(indentations, ' ') << "<NODE name=\"" << escapeXML(trace.name) << "\" description=\"" << escapeXML(d) << "\">\n";
          indentations += 2;
        }
        else
        {
          indentations -= 2;
          os << std::string(indentations, ' ') << "</NODE>\n";
        }
      }

      if (param_it->value.valueType() != ParamValue::EMPTY_VALUE)
      {
        // we create a temporary copy of the tag list, since we remove certain tags while writing,
        // that will be represented differently in the xml
        std::set<std::string> tag_list = param_it->tags;
        ParamValue::ValueType value_type = param_it->value.valueType();
        bool stringParamIsFlag = false;

        if (value_type < ParamValue::STRING_LIST)
        {
          os << std::string(indentations, ' ') << "<ITEM name=\"" << escapeXML(param_it->name) << R"(" value=")";
        }
        else
        {
          os << std::string(indentations, ' ') << "<ITEMLIST name=\"" << escapeXML(param_it->name);
        }

        switch (value_type)
        {
          case ParamValue::INT_VALUE:
            os << param_it->value.toString() << R"(" type="int")";
            break;
          case ParamValue::DOUBLE_VALUE:
            os << param_it->value.toString() << R"(" type="double")";
            break;
          case ParamValue::STRING_VALUE:
            if (tag_list.find(TOPPBase::TAG_INPUT_FILE) != tag_list.end())
            {
              os << escapeXML(param_it->value.toString()) << R"(" type="input-file")";
              tag_list.erase(TOPPBase::TAG_INPUT_FILE);
            }
            else if (tag_list.find(TOPPBase::TAG_OUTPUT_FILE) != tag_list.end())
            {
              os << escapeXML(param_it->value.toString()) << R"(" type="output-file")";
              tag_list.erase(TOPPBase::TAG_OUTPUT_FILE);
            }
            else if (tag_list.find(TOPPBase::TAG_OUTPUT_DIR) != tag_list.end())
            {
              os << escapeXML(param_it->value.toString()) << R"(" type="output-dir")";
              tag_list.erase(TOPPBase::TAG_OUTPUT_DIR);
            }
            else if (tag_list.find(TOPPBase::TAG_OUTPUT_PREFIX) != tag_list.end())
            {
              os << escapeXML(param_it->value.toString()) << R"(" type="output-prefix")";
              tag_list.erase(TOPPBase::TAG_OUTPUT_PREFIX);
            }
            else if (param_it->valid_strings.size() == 2 && param_it->valid_strings[0] == "true" && param_it->valid_strings[1] == "false" && param_it->value == "false")
            {
              stringParamIsFlag = true;
              os << param_it->value.toString() << R"(" type="bool")";
            }
            else
            {
              std::string value = param_it->value.toString();
              if (value.find('\t') != std::string::npos)
              {
                replace(value, '\t', "&#x9;");
              }
              os << escapeXML(value) << R"(" type="string")";
            }
            break;
          case ParamValue::STRING_LIST:
            if (tag_list.find(TOPPBase::TAG_INPUT_FILE) != tag_list.end())
            {
              os << R"(" type="input-file")";
              tag_list.erase(TOPPBase::TAG_INPUT_FILE);
            }
            else if (tag_list.find(TOPPBase::TAG_OUTPUT_FILE) != tag_list.end())
            {
              os << R"(" type="output-file")";
              tag_list.erase(TOPPBase::TAG_OUTPUT_FILE);
            }
            else
            {
              os << R"(" type="string")";
            }
            break;
          case ParamValue::INT_LIST:
            os << R"(" type="int")";
            break;
          case ParamValue::DOUBLE_LIST:
            os << R"(" type="double")";
            break;
          default:
            break;
        }

        std::string description = param_it->description;
        replace(description, '\n', "#br#");

        os << " description=\"" << escapeXML(description) << "\"";

        if (tag_list.find("required") != tag_list.end())
        {
          os << R"( required="true")";
          tag_list.erase("required");
        }
        else
        {
          os << R"( required="false")";
        }

        if (tag_list.find("advanced") != tag_list.end())
        {
          os << R"( advanced="true")";
          tag_list.erase("advanced");
        }
        else
        {
          os << R"( advanced="false")";
        }

        if (!tag_list.empty())
        {
          std::string list;
          for (auto& tag : tag_list)
          {
            if (!list.empty())
              list += ",";
            list += tag;
          }
          os << " tags=\"" << escapeXML(list) << "\"";
        }

        if (!stringParamIsFlag)
        {
          std::string restrictions;
          switch (value_type)
          {
            case ParamValue::INT_VALUE:
            case ParamValue::INT_LIST: {
              bool min_set = (param_it->min_int != -std::numeric_limits<int>::max());
              bool max_set = (param_it->max_int != std::numeric_limits<int>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions += std::to_string(param_it->min_int);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions += std::to_string(param_it->max_int);
                }
              }
            }
            break;

            case ParamValue::DOUBLE_VALUE:
            case ParamValue::DOUBLE_LIST: {
              bool min_set = (param_it->min_float != -std::numeric_limits<double>::max());
              bool max_set = (param_it->max_float != std::numeric_limits<double>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions += std::to_string(param_it->min_float);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions += std::to_string(param_it->max_float);
                }
              }
            }
            break;
            case ParamValue::STRING_VALUE:
            case ParamValue::STRING_LIST:
              if (!param_it->valid_strings.empty())
              {
                restrictions = param_it->valid_strings.front();
                for (auto it = param_it->valid_strings.begin() + 1, end = param_it->valid_strings.end(); it != end; ++it)
                {
                  restrictions += ",";
                  restrictions += *it;
                }
              }
              break;
            default:
              break;
          }

          if (!restrictions.empty())
          {
            if (param_it->tags.find("input file") != param_it->tags.end() || param_it->tags.find("output file") != param_it->tags.end() || param_it->tags.find("output prefix") != param_it->tags.end())
            {
              os << " supported_formats=\"" << escapeXML(restrictions) << "\"";
            }
            else
            {
              os << " restrictions=\"" << escapeXML(restrictions) << "\"";
            }
          }
        }

        if (value_type < ParamValue::STRING_LIST)
        {
          os << " />\n";
        }
        else
        {
          os << " >\n";
        }

        switch (value_type)
        {
          case ParamValue::STRING_LIST:
            for (auto item : static_cast<const std::vector<std::string>&>(param_it->value))
            {
              if (item.find('\t') != std::string::npos)
              {
                replace(item, '\t', "&#x9;");
              }
              os << std::string(indentations + 2, ' ') << "<LISTITEM value=\"" << escapeXML(item) << "\"/>\n";
            }
            break;
          case ParamValue::INT_LIST:
            for (int item : static_cast<const std::vector<int>&>(param_it->value))
            {
              os << std::string(indentations + 2, ' ') << "<LISTITEM value=\"" << item << "\"/>\n";
            }
            break;
          case ParamValue::DOUBLE_LIST:
            for (double item : static_cast<const std::vector<double>&>(param_it->value))
            {
              os << std::string(indentations + 2, ' ') << "<LISTITEM value=\"" << item << "\"/>\n";
            }
            break;
          default:
            break;
        }

        if (value_type > ParamValue::DOUBLE_VALUE && value_type)
        {
          os << std::string(indentations, ' ') << "</ITEMLIST>\n";
        }
      }
    }

    if (param.begin() != param.end())
    {
      for ([[maybe_unused]] auto& trace : param_it.getTrace())
      {
        indentations -= 2;
        os << std::string(indentations, ' ') << "</NODE>\n";
      }
    }

    os << "</PARAMETERS>\n";
    os << "</tool>" << std::endl; // forces a flush
  }

  std::string ParamCTDFile::escapeXML(const std::string& to_escape)
  {
    std::string copy = to_escape;
    if (copy.find('&') != std::string::npos)
      replace(copy, '&', "&amp;");
    if (copy.find('>') != std::string::npos)
      replace(copy, '>', "&gt;");
    if (copy.find('"') != std::string::npos)
      replace(copy, '"', "&quot;");
    if (copy.find('<') != std::string::npos)
      replace(copy, '<', "&lt;");
    if (copy.find('\'') != std::string::npos)
      replace(copy, '\'', "&apos;");

    return copy;
  }

  void ParamCTDFile::replace(std::string& replace_in, char to_replace, const std::string& replace_with)
  {
    for (size_t i = 0; i < replace_in.size(); ++i)
    {
      if (replace_in[i] == to_replace)
      {
        replace_in = replace_in.substr(0, i) + replace_with + replace_in.substr(i + 1);
        i += replace_with.size();
      }
    }
  }
} // namespace OpenMS
