// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <iostream>
#include <fstream>
#include <algorithm>

namespace OpenMS
{

  String writeXMLEscape(const String& to_escape)
  {
    return Internal::XMLHandler::writeXMLEscape(to_escape);
  }

  ParamXMLFile::ParamXMLFile() :
    XMLFile("/SCHEMAS/Param_1_8_0.xsd", "1.8.0")
  {
  }

  void ParamXMLFile::store(const String& filename, const Param& param) const
  {
    //open file
    std::ofstream os_;
    std::ostream* os_ptr;
    if (filename != "-")
    {
      os_.open(filename.c_str(), std::ofstream::out);
      if (!os_)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      os_ptr = &os_;
    }
    else
    {
      os_ptr = &std::cout;
    }

    //write to file stream
    writeXMLToStream(os_ptr, param);

    os_.close();
  }

  void ParamXMLFile::writeXMLToStream(std::ostream* os_ptr, const Param& param) const
  {
    // Note: For a long time the handling of 'getTrace()' was vulnerable to an unpruned tree (a path of nodes, but no entries in them), i.e.
    //       too many closing tags are written to the INI file, but no opening ones.
    //       This never mattered here, as removeAll() was fixed to prune the tree.
    // TODO: Nowadays this should be fixed and removeAll() might not be necessary.

    std::ostream& os = *os_ptr;

    os.precision(writtenDigits<double>(0.0));

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<PARAMETERS version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/Param_1_8_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    String indentation = "  ";
    Param::ParamIterator it = param.begin();
    while (it != param.end())
    {
      //write opened/closed nodes
      const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (const Param::ParamIterator::TraceInfo& it2 : trace)
      {
        if (it2.opened) //opened node
        {
          String d = it2.description;
          //d.substitute('"','\'');
          d.substitute("\n", "#br#");
          //d.substitute("<","&lt;");
          //d.substitute(">","&gt;");
          os << indentation  << "<NODE name=\"" << writeXMLEscape(it2.name) << "\" description=\"" << writeXMLEscape(d) << "\">" << "\n";
          indentation += "  ";
        }
        else //closed node
        {
          indentation.resize(indentation.size() - 2);
          os << indentation << "</NODE>" << "\n";
        }
      }

      //write item
      if (it->value.valueType() != ParamValue::EMPTY_VALUE)
      {
        // we create a temporary copy of the tag list, since we remove certain tags while writing,
        // that will be represented differently in the xml
        std::set<std::string> tag_list = it->tags;
        ParamValue::ValueType value_type = it->value.valueType();
        bool stringParamIsFlag = false;

        //write opening tag
        switch (value_type)
        {
        case ParamValue::INT_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << R"(" type="int")";
          break;

        case ParamValue::DOUBLE_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << R"(" type="double")";
          break;

        case ParamValue::STRING_VALUE:
          if (tag_list.find("input file") != tag_list.end())
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="input-file")";
            tag_list.erase("input file");
          }
          else if (tag_list.find("output file") != tag_list.end())
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="output-file")";
            tag_list.erase("output file");
          }
          else if (tag_list.find("output prefix") != tag_list.end())
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << R"(" type="output-prefix")";
            tag_list.erase("output prefix");
          }

          else if (it->valid_strings.size() == 2 &&
          it->valid_strings[0] == "true" && it->valid_strings[1] == "false" &&
          it->value == "false")
          {
            stringParamIsFlag = true;
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << R"(" type="bool")";
          }
          else
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << R"(" type="string")";
          }
          break;

        case ParamValue::STRING_LIST:
          if (tag_list.find("input file") != tag_list.end())
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="input-file")";
            tag_list.erase("input file");
          }
          else if (tag_list.find("output file") != tag_list.end())
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="output-file")";
            tag_list.erase("output file");
          }
          else
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="string")";
          }
          break;

        case ParamValue::INT_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="int")";
          break;

        case ParamValue::DOUBLE_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << R"(" type="double")";
          break;

        default:
          break;
        }

        //replace all critical characters in description
        String d = it->description;
        //d.substitute("\"","'");
        d.substitute("\n", "#br#");
        //d.substitute("<","&lt;");
        //d.substitute(">","&gt;");
        os << " description=\"" << writeXMLEscape(d) << "\"";

        // required
        if (tag_list.find("required") != tag_list.end())
        {
          os << " required=\"true\"";
          tag_list.erase("required");
        }
        else
        {
          os << " required=\"false\"";
        }

        // advanced
        if (tag_list.find("advanced") != tag_list.end())
        {
          os << " advanced=\"true\"";
          tag_list.erase("advanced");
        }
        else
        {
          os << " advanced=\"false\"";
        }

        // tags
        if (!tag_list.empty())
        {
          String list;
          for (std::set<std::string>::const_iterator tag_it = tag_list.begin(); tag_it != tag_list.end(); ++tag_it)
          {
            if (!list.empty())
              list += ",";
            list += *tag_it;
          }
          os << " tags=\"" << writeXMLEscape(list) << "\"";
        }

        //restrictions
        // for boolean Flags they are implicitly given
        if (!stringParamIsFlag)
        {
          String restrictions = "";
          switch (value_type)
          {
            case ParamValue::INT_VALUE:
            case ParamValue::INT_LIST:
            {
              bool min_set = (it->min_int != -std::numeric_limits<Int>::max());
              bool max_set = (it->max_int != std::numeric_limits<Int>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions += String(it->min_int);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions += String(it->max_int);
                }
              }
            }
              break;

            case ParamValue::DOUBLE_VALUE:
            case ParamValue::DOUBLE_LIST:
            {
              bool min_set = (it->min_float != -std::numeric_limits<double>::max());
              bool max_set = (it->max_float != std::numeric_limits<double>::max());
              if (max_set || min_set)
              {
                if (min_set)
                {
                  restrictions += String(it->min_float);
                }
                restrictions += ':';
                if (max_set)
                {
                  restrictions += String(it->max_float);
                }
              }
            }
              break;

            case ParamValue::STRING_VALUE:
            case ParamValue::STRING_LIST:
              if (!it->valid_strings.empty())
              {
                restrictions.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ",");
              }
              break;

            default:
              break;
          }
          // for files we store the restrictions as supported_formats
          if (!restrictions.empty())
          {
            if (it->tags.find("input file") != it->tags.end() 
              || it->tags.find("output file") != it->tags.end()
              || it->tags.find("output prefix") != it->tags.end())
            {
              os << " supported_formats=\"" << writeXMLEscape(restrictions) << "\"";
            }
            else
            {
              os << " restrictions=\"" << writeXMLEscape(restrictions) << "\"";
            }
          }
        }

        //finish opening tag
        switch (value_type)
        {
        case ParamValue::INT_VALUE:
        case ParamValue::DOUBLE_VALUE:
        case ParamValue::STRING_VALUE:
          os << " />" <<  "\n";
          break;

        case ParamValue::STRING_LIST:
        {
          os << ">" <<  "\n";
          const std::vector<std::string>& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << Internal::encodeTab(writeXMLEscape(list[i])) << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        case ParamValue::INT_LIST:
        {
          os << ">" <<  "\n";
          const IntList& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        case ParamValue::DOUBLE_LIST:
        {
          os << ">" <<  "\n";
          const DoubleList& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << list[i] << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        default:
          break;
        }
      }
      ++it;
    }

    // if we had tags ...
    if (param.begin() != param.end())
    {
      //close remaining tags
      const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        Size ss = indentation.size();
        indentation.resize(ss - 2);
        os << indentation << "</NODE>" << "\n";
      }
    }

    os << "</PARAMETERS>" << std::endl; // forces a flush
  }

  void ParamXMLFile::load(const String& filename, Param& param)
  {
    Internal::ParamXMLHandler handler(param, filename, schema_version_);
    parse_(filename, &handler);
  }

}
