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
    XMLFile("/SCHEMAS/Param_1_7_0.xsd", "1.7.0")
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
    // hint: the handling of 'getTrace()' is vulnerable to an unpruned tree (a path of nodes, but no entries in them), i.e.
    //       too many closing tags are written to the INI file, but no opening ones.
    //       This currently cannot happen, as removeAll() was fixed to prune the tree, just keep it in mind.

    std::ostream& os = *os_ptr;

    os.precision(writtenDigits<double>(0.0));

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<PARAMETERS version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/Param_1_7_0.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    String indentation = "  ";
    Param::ParamIterator it = param.begin();
    while (it != param.end())
    {
      //write opened/closed nodes
      const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        if (it2->opened) //opened node
        {
          String d = it2->description;
          //d.substitute('"','\'');
          d.substitute("\n", "#br#");
          //d.substitute("<","&lt;");
          //d.substitute(">","&gt;");
          os << indentation  << "<NODE name=\"" << writeXMLEscape(it2->name) << "\" description=\"" << writeXMLEscape(d) << "\">" << "\n";
          indentation += "  ";
        }
        else //closed node
        {
          indentation.resize(indentation.size() - 2);
          os << indentation << "</NODE>" << "\n";
        }
      }

      //write item
      if (it->value.valueType() != DataValue::EMPTY_VALUE)
      {
        // we create a temporary copy of the tag list, since we remove certain tags while writing,
        // that will be represented differently in the xml
        std::set<String> tag_list = it->tags;
        DataValue::DataType value_type = it->value.valueType();
        bool stringParamIsFlag = false;

        //write opening tag
        switch (value_type)
        {
        case DataValue::INT_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"int\"";
          break;

        case DataValue::DOUBLE_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"double\"";
          break;

        case DataValue::STRING_VALUE:
          if (tag_list.find("input file") != tag_list.end())
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << "\" type=\"input-file\"";
            tag_list.erase("input file");
          }
          else if (tag_list.find("output file") != tag_list.end())
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << "\" type=\"output-file\"";
            tag_list.erase("output file");
          }
          else if (it->valid_strings.size() == 2 &&
          it->valid_strings[0] == "true" && it->valid_strings[1] == "false" &&
          it->value == "false")
          {
            stringParamIsFlag = true;
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << "\" type=\"bool\"";
          }
          else
          {
            os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << Internal::encodeTab(writeXMLEscape(it->value.toString())) << "\" type=\"string\"";
          }
          break;

        case DataValue::STRING_LIST:
          if (tag_list.find("input file") != tag_list.end())
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"input-file\"";
            tag_list.erase("input file");
          }
          else if (tag_list.find("output file") != tag_list.end())
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"output-file\"";
            tag_list.erase("output file");
          }
          else
          {
            os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"string\"";
          }
          break;

        case DataValue::INT_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"int\"";
          break;

        case DataValue::DOUBLE_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"double\"";
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
          for (std::set<String>::const_iterator tag_it = tag_list.begin(); tag_it != tag_list.end(); ++tag_it)
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
            case DataValue::INT_VALUE:
            case DataValue::INT_LIST:
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

            case DataValue::DOUBLE_VALUE:
            case DataValue::DOUBLE_LIST:
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

            case DataValue::STRING_VALUE:
            case DataValue::STRING_LIST:
              if (it->valid_strings.size() != 0)
              {
                restrictions.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ",");
              }
              break;

            default:
              break;
          }
          // for files we store the restrictions as supported_formats
          if (restrictions != "")
          {
            if (it->tags.find("input file") != it->tags.end() || it->tags.find("output file") != it->tags.end())
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
        case DataValue::INT_VALUE:
        case DataValue::DOUBLE_VALUE:
        case DataValue::STRING_VALUE:
          os << " />" <<  "\n";
          break;

        case DataValue::STRING_LIST:
        {
          os << ">" <<  "\n";
          const StringList& list = it->value;
          for (Size i = 0; i < list.size(); ++i)
          {
            os << indentation << "  <LISTITEM value=\"" << Internal::encodeTab(writeXMLEscape(list[i])) << "\"/>" << "\n";
          }
          os << indentation << "</ITEMLIST>" << "\n";
        }
        break;

        case DataValue::INT_LIST:
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

        case DataValue::DOUBLE_LIST:
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
