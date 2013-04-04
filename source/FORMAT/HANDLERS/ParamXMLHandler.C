// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <iostream>

#include <xercesc/sax2/Attributes.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    ParamXMLHandler::ParamXMLHandler(Param& param, const String& filename, const String& version) :
      XMLHandler(filename, version),
      param_(param)
    {
    }

    ParamXMLHandler::~ParamXMLHandler()
    {
    }

    void ParamXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
    {
      static const XMLCh* s_restrictions = xercesc::XMLString::transcode("restrictions");
      static const XMLCh* s_supported_formats = xercesc::XMLString::transcode("supported_formats");

      String element = sm_.convert(qname);
      if (element == "ITEM")
      {
        //parse value/type
        String type = attributeAsString_(attributes, "type");
        String name = path_ + attributeAsString_(attributes, "name");
        String value = attributeAsString_(attributes, "value");
        //parse description, if present
        String description;
        optionalAttributeAsString_(description, attributes, "description");
        description.substitute("#br#", "\n");
        //tags
        String tags_string;
        optionalAttributeAsString_(tags_string, attributes, "tags");
        StringList tags = StringList::create(tags_string);
        //advanced (for downward compatibility with old Param files)
        String advanced_string;
        optionalAttributeAsString_(advanced_string, attributes, "advanced");
        if (advanced_string == "true")
        {
          tags.push_back("advanced");
        }
        //type
        if (type == "int")
        {
          param_.setValue(name, asInt_(value), description, tags);
        }
        else if (type == "string")
        {
          param_.setValue(name, value, description, tags);
        }
        else if (type == "float" || type == "double")
        {
          param_.setValue(name, asDouble_(value), description, tags);
        }
        else
        {
          warning(LOAD, String("Ignoring entry '") + name + "' because of unknown type '" + type + "'");
        }

        //restrictions
        Int restrictions_index = attributes.getIndex(s_restrictions);
        if (restrictions_index != -1)
        {
          String value = sm_.convert(attributes.getValue(restrictions_index));
          std::vector<String> parts;
          if (type == "int")
          {
            value.split(':', parts);
            if (parts.size() != 2)
              value.split('-', parts);  //for downward compatibility
            if (parts.size() == 2)
            {
              if (parts[0] != "")
              {
                param_.setMinInt(name, parts[0].toInt());
              }
              if (parts[1] != "")
              {
                param_.setMaxInt(name, parts[1].toInt());
              }
            }
            else
            {
              warning(LOAD, "ITEM " + name + " has an empty restrictions attribute.");
            }
          }
          else if (type == "string")
          {
            value.split(',', parts);
            param_.setValidStrings(name, parts);
          }
          else if (type == "float" || type == "double")
          {
            value.split(':', parts);
            if (parts.size() != 2)
              value.split('-', parts);  //for downward compatibility
            if (parts.size() == 2)
            {
              if (parts[0] != "")
              {
                param_.setMinFloat(name, parts[0].toDouble());
              }
              if (parts[1] != "")
              {
                param_.setMaxFloat(name, parts[1].toDouble());
              }
            }
            else
            {
              warning(LOAD, "ITEM " + name + " has an empty restrictions attribute.");
            }
          }
        }

        // check for supported_formats -> supported_formats overwrites restrictions in case of files
        if ((tags.contains("input file") || tags.contains("output file")) && type == "string")
        {
          Int supported_formats_index = attributes.getIndex(s_supported_formats);
          if (supported_formats_index != -1)
          {
            String value = sm_.convert(attributes.getValue(supported_formats_index));
            std::vector<String> parts;

            value.split(',', parts);
            param_.setValidStrings(name, parts);
          }
        }

      }
      else if (element == "NODE")
      {
        //parse name
        String name = attributeAsString_(attributes, "name");
        open_tags_.push_back(name);
        path_ += name + ":";
        //parse description
        String description;
        optionalAttributeAsString_(description, attributes, "description");
        if (description != "")
        {
          description.substitute("#br#", "\n");
          descriptions_[path_.chop(1)] = description;
        }
      }
      else if (element == "ITEMLIST")
      {
        //parse name/type
        list_.type = attributeAsString_(attributes, "type");
        list_.name = path_ + attributeAsString_(attributes, "name");
        //parse description, if present
        list_.description = "";
        optionalAttributeAsString_(list_.description, attributes, "description");
        list_.description.substitute("#br#", "\n");
        //tags
        String tags_string;
        optionalAttributeAsString_(tags_string, attributes, "tags");
        list_.tags = StringList::create(tags_string);
        //advanced (for downward compatibility with old Param files)
        String advanced_string;
        optionalAttributeAsString_(advanced_string, attributes, "advanced");
        if (advanced_string == "true")
        {
          list_.tags.push_back("advanced");
        }
        list_.restrictions_index = attributes.getIndex(s_restrictions);
        if (list_.restrictions_index != -1)
        {
          list_.restrictions = sm_.convert(attributes.getValue(list_.restrictions_index));
        }

        // check for supported_formats -> supported_formats overwrites restrictions in case of files
        if ((list_.tags.contains("input file") || list_.tags.contains("output file")) && list_.type == "string")
        {
          Int supported_formats_index = attributes.getIndex(s_supported_formats);
          if (supported_formats_index != -1)
          {
            list_.restrictions_index = supported_formats_index;
            list_.restrictions = sm_.convert(attributes.getValue(list_.restrictions_index));
          }
        }
      }
      else if (element == "LISTITEM")
      {
        if (list_.type == "string")
        {
          list_.stringlist.push_back(attributeAsString_(attributes, "value"));
        }
        else if (list_.type == "int")
        {
          list_.intlist.push_back(asInt_(attributeAsString_(attributes, "value")));
        }
        else if (list_.type == "float" || list_.type == "double")
        {
          list_.doublelist.push_back(asDouble_(attributeAsString_(attributes, "value")));
        }

      }
      else if (element == "PARAMETERS")
      {
        //check file version against schema version
        String file_version = "";
        optionalAttributeAsString_(file_version, attributes, "version");
        if (file_version == "")
          file_version = "1.0";                       //default version is 1.0
        if (file_version.toDouble() > version_.toDouble())
        {
          warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
        }
      }
    }

    void ParamXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      String element = sm_.convert(qname);
      if (element == "NODE")
      {
        open_tags_.pop_back();
        //renew path
        path_ = "";
        for (vector<String>::iterator it = open_tags_.begin(); it != open_tags_.end(); ++it)
        {
          path_ += *it + ":";
        }
      }
      else if (element == "ITEMLIST")
      {
        std::vector<String> parts;

        if (list_.type == "string")
        {
          param_.setValue(list_.name, list_.stringlist, list_.description, list_.tags);
          if (list_.restrictions_index != -1)
          {
            list_.restrictions.split(',', parts);
            param_.setValidStrings(list_.name, parts);
          }
        }
        else if (list_.type == "int")
        {
          param_.setValue(list_.name, list_.intlist, list_.description, list_.tags);
          if (list_.restrictions_index != -1)
          {
            list_.restrictions.split(':', parts);
            if (parts.size() != 2)
              list_.restrictions.split('-', parts);  //for downward compatibility
            if (parts.size() == 2)
            {
              if (parts[0] != "")
              {
                param_.setMinInt(list_.name, parts[0].toInt());
              }
              if (parts[1] != "")
              {
                param_.setMaxInt(list_.name, parts[1].toInt());
              }
            }
            else
            {
              warning(LOAD, "ITEMLIST " + list_.name + " has an empty restrictions attribute.");
            }
          }
        }
        else if (list_.type == "float" || list_.type == "double")
        {
          param_.setValue(list_.name, list_.doublelist, list_.description, list_.tags);
          if (list_.restrictions_index != -1)
          {
            list_.restrictions.split(':', parts);
            if (parts.size() != 2)
              list_.restrictions.split('-', parts);  //for downward compatibility
            if (parts.size() == 2)
            {
              if (parts[0] != "")
              {
                param_.setMinFloat(list_.name, parts[0].toDouble());
              }
              if (parts[1] != "")
              {
                param_.setMaxFloat(list_.name, parts[1].toDouble());
              }
            }
            else
            {
              warning(LOAD, "ITEMLIST " + list_.name + " has an empty restrictions attribute.");
            }
          }
        }
        else
        {
          warning(LOAD, String("Ignoring list entry '") + list_.name + "' because of unknown type '" + list_.type + "'");
        }
        list_.stringlist.clear();
        list_.intlist.clear();
        list_.doublelist.clear();
      }
      else if (element == "PARAMETERS")
      {
        //set all descriptions (now the nodes exist...)
        for (map<String, String>::const_iterator it = descriptions_.begin(); it != descriptions_.end(); ++it)
        {
          param_.setSectionDescription(it->first, it->second);
        }
        descriptions_.clear();
      }
    }

  } // namespace Internal
} // namespace OpenMS
