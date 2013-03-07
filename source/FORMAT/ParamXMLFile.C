// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <fstream>
#include <set>

using namespace OpenMS::Internal;
using namespace std;

namespace OpenMS
{
  ParamXMLFile::ParamXMLFile() :
    XMLFile("/SCHEMAS/Param_1_4.xsd", "1.4")
  {
  }

  void ParamXMLFile::store(const String& filename, const Param& param) const
  {
    //open file
    ofstream os_;
    ostream* os_ptr;
    if (filename != "-")
    {
      os_.open(filename.c_str(), std::ofstream::out);
      if (!os_)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
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

  void ParamXMLFile::writeXMLToStream(ostream* os_ptr, const Param& param) const
  {
    // hint: the handling of 'getTrace()' is vulnerable to an unpruned tree (a path of nodes, but no entries in them), i.e.
    //       too many closing tags are written to the INI file, but no openening ones.
    //       This currently cannot happen, as removeAll() was fixed to prune the tree, just keep it in mind.

    ostream& os = *os_ptr;

    os.precision(writtenDigits<DoubleReal>());

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<PARAMETERS version=\"" << getVersion() << "\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_4.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    String indentation = "  ";
    Param::ParamIterator it = param.begin();
    while (it != param.end())
    {
      //write opened/closed nodes
      const vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
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
        else         //closed node
        {
          indentation.resize(indentation.size() - 2);
          os << indentation << "</NODE>" << "\n";
        }
      }

      //write item
      if (it->value.valueType() != DataValue::EMPTY_VALUE)
      {
        DataValue::DataType value_type = it->value.valueType();
        //write opening tag
        switch (value_type)
        {
        case DataValue::INT_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"int\"";
          break;

        case DataValue::DOUBLE_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << it->value.toString() << "\" type=\"float\"";
          break;

        case DataValue::STRING_VALUE:
          os << indentation << "<ITEM name=\"" << writeXMLEscape(it->name) << "\" value=\"" << writeXMLEscape(it->value.toString()) << "\" type=\"string\"";
          break;

        case DataValue::STRING_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"string\"";
          break;

        case DataValue::INT_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"int\"";
          break;

        case DataValue::DOUBLE_LIST:
          os << indentation << "<ITEMLIST name=\"" << writeXMLEscape(it->name) << "\" type=\"float\"";
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

        //tags
        if (!it->tags.empty())
        {
          String list;
          for (set<String>::const_iterator tag_it = it->tags.begin(); tag_it != it->tags.end(); ++tag_it)
          {
            if (!list.empty())
              list += ",";
            list += *tag_it;
          }
          os << " tags=\"" << writeXMLEscape(list) << "\"";
        }

        //restrictions
        String restrictions = "";
        switch (value_type)
        {
        case DataValue::INT_VALUE:
        case DataValue::INT_LIST:
        {
          bool min_set = (it->min_int != -numeric_limits<Int>::max());
          bool max_set = (it->max_int != numeric_limits<Int>::max());
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
          bool min_set = (it->min_float != -numeric_limits<DoubleReal>::max());
          bool max_set = (it->max_float != numeric_limits<DoubleReal>::max());
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
            os << indentation << "  <LISTITEM value=\"" << writeXMLEscape(list[i]) << "\"/>" << "\n";
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
      const vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
      for (vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
      {
        Size ss = indentation.size();
        indentation.resize(ss - 2);
        os << indentation << "</NODE>" << "\n";
      }
    }

    os << "</PARAMETERS>" << endl; // forces a flush
  }

  void ParamXMLFile::load(const String& filename, Param& param)
  {
    ParamXMLHandler handler(param, filename, schema_version_);
    parse_(filename, &handler);
  }

}
