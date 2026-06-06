// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

using namespace xercesc;
using namespace std;

namespace OpenMS::Internal
{


    ParamXMLHandler::ParamXMLHandler(Param& param, const std::string& filename, const std::string& version) :
      XMLHandler(filename, version),
      param_(param)
    {
    }

    ParamXMLHandler::~ParamXMLHandler()
    = default;

    void ParamXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
    {
      static const XMLCh* s_restrictions = xercesc::XMLString::transcode("restrictions");
      static const XMLCh* s_supported_formats = xercesc::XMLString::transcode("supported_formats");

      std::string element = sm_.convert(qname);
      if (element == "ITEM")
      {
        //parse value/type
        std::string type = attributeAsString_(attributes, "type");
        std::string name = path_ + attributeAsString_(attributes, "name");
        std::string value = attributeAsString_(attributes, "value");

        //parse description, if present
        std::string description;
        optionalAttributeAsString_(description, attributes, "description");
        StringUtils::substitute(description, "#br#", "\n");

        //tags
        std::string tags_string;
        optionalAttributeAsString_(tags_string, attributes, "tags");
        std::vector<std::string> tags = ListUtils::create<std::string>(tags_string);

        //advanced
        std::string advanced_string;
        optionalAttributeAsString_(advanced_string, attributes, "advanced");
        if (advanced_string == "true")
        {
          tags.emplace_back("advanced");
        }

        //required
        std::string required_string;
        optionalAttributeAsString_(advanced_string, attributes, "required");
        if (advanced_string == "true")
        {
          tags.emplace_back("required");
        }

        //type
        if (type == "int")
        {
          param_.setValue(name, asInt_(value), description, tags);
        }
        // since 1.7 it supports a separate bool type
        else if (type == "bool" || type == "string")
        {
          param_.setValue(name, value, description, tags);
        }
        // since param v1.6.2 we support explicitly naming input/output files as types
        else if (type == "input-file")
        {
          tags.emplace_back("input file");
          param_.setValue(name, value, description, tags);
        }
        else if (type == "output-file")
        {
          tags.emplace_back("output file");          
          param_.setValue(name, value, description, tags);
        }
        else if (type == "output-prefix")
        {
          tags.emplace_back("output prefix");          
          param_.setValue(name, value, description, tags);
        }
        else if (type == "float" || type == "double")
        {
          param_.setValue(name, asDouble_(value), description, tags);
        }
        else
        {
          warning(LOAD,StringUtils::toStr("Ignoring entry '") + name + "' because of unknown type '" + type + "'");
        }

        //restrictions

        // we internally handle bool parameters as strings with restrictions true/false
        if (type == "bool")
        {
          param_.setValidStrings(name, {"true","false"});
        }
        else
        {
          // parse restrictions if present
          Int restrictions_index = attributes.getIndex(s_restrictions);
          if (restrictions_index != -1)
          {
            std::string val = sm_.convert(attributes.getValue(restrictions_index));
            std::vector<std::string> parts;
            if (type == "int")
            {
              StringUtils::split(val, ':', parts);
              if (parts.size() != 2)
                StringUtils::split(val, '-', parts); //for downward compatibility
              if (parts.size() == 2)
              {
                if (!parts[0].empty())
                {
                  param_.setMinInt(name, StringUtils::toInt32(parts[0]));
                }
                if (!parts[1].empty())
                {
                  param_.setMaxInt(name, StringUtils::toInt32(parts[1]));
                }
              }
              else
              {
                warning(LOAD, "ITEM " + name + " has an empty restrictions attribute.");
              }
            }
            else if (type == "string")
            {
              StringUtils::split(val, ', ', parts);
              param_.setValidStrings(name, ListUtils::create<std::string>(parts));
            }
            else if (type == "float" || type == "double")
            {
              StringUtils::split(val, ':', parts);
              if (parts.size() != 2)
              {
                StringUtils::split(val, '-', parts); //for downward compatibility
              }
              if (parts.size() == 2)
              {
                if (!parts[0].empty())
                {
                  param_.setMinFloat(name, StringUtils::toDouble(parts[0]));
                }
                if (!parts[1].empty())
                {
                  param_.setMaxFloat(name, StringUtils::toDouble(parts[1]));
                }
              }
              else
              {
                warning(LOAD, "ITEM " + name + " has an empty restrictions attribute.");
              }
            }
          }
        }


        // check for supported_formats -> supported_formats overwrites restrictions in case of files
        if ((ListUtils::contains(tags, "input file") || ListUtils::contains(tags, "output file")) && (type == "string" || type == "input-file" || type == "output-file"))
        {
          Int supported_formats_index = attributes.getIndex(s_supported_formats);
          if (supported_formats_index != -1)
          {
            std::string val = sm_.convert(attributes.getValue(supported_formats_index));
            std::vector<std::string> parts;

            StringUtils::split(val, ', ', parts);
            param_.setValidStrings(name, ListUtils::create<std::string>(parts));
          }
        }

      }
      else if (element == "NODE")
      {
        //parse name
        std::string name = attributeAsString_(attributes, "name");
        open_tags_.push_back(name);
        path_ += name + ":";
        //parse description
        std::string description;
        optionalAttributeAsString_(description, attributes, "description");
        if (!description.empty())
        {
          StringUtils::substitute(description, "#br#", "\n");
        }
        param_.addSection(StringUtils::chop(path_, 1), description);
      }
      else if (element == "ITEMLIST")
      {
        //tags
        std::string tags_string;
        optionalAttributeAsString_(tags_string, attributes, "tags");
        list_.tags = ListUtils::create<std::string>(tags_string);
                
        //parse name/type
        list_.type = attributeAsString_(attributes, "type");
        // handle in-/output file correctly
        if (list_.type == "input-file")
        {
          list_.type = "string";
          list_.tags.emplace_back("input file");
        }
        else if (list_.type == "output-file")
        {
          list_.type = "string";
          list_.tags.emplace_back("output file");
        }

        list_.name = path_ + attributeAsString_(attributes, "name");
        
        //parse description, if present
        list_.description = "";
        optionalAttributeAsString_(list_.description, attributes, "description");
        StringUtils::substitute(list_.description, "#br#", "\n");
        
        //advanced
        std::string advanced_string;
        optionalAttributeAsString_(advanced_string, attributes, "advanced");
        if (advanced_string == "true")
        {
          list_.tags.emplace_back("advanced");
        }

        //advanced
        std::string required_string;
        optionalAttributeAsString_(required_string, attributes, "required");
        if (required_string == "true")
        {
          list_.tags.emplace_back("required");
        }
        
        list_.restrictions_index = attributes.getIndex(s_restrictions);
        if (list_.restrictions_index != -1)
        {
          list_.restrictions = sm_.convert(attributes.getValue(list_.restrictions_index));
        }

        // check for supported_formats -> supported_formats overwrites restrictions in case of files
        if ((ListUtils::contains(list_.tags, "input file") || ListUtils::contains(list_.tags, "output file")) && list_.type == "string")
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
        std::string file_version = "";
        optionalAttributeAsString_(file_version, attributes, "version");

        // default version is 1.0
        if (file_version.empty())
        {
          file_version = "1.0";
        }
        VersionInfo::VersionDetails file_version_details = VersionInfo::VersionDetails::create(file_version);
        VersionInfo::VersionDetails parser_version = VersionInfo::VersionDetails::create(version_);

        if (file_version_details > parser_version)
        {
          warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
        }
      }
    }

    void ParamXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      std::string element = sm_.convert(qname);
      if (element == "NODE")
      {
        open_tags_.pop_back();
        //renew path
        path_ = "";
        for (vector<std::string>::iterator it = open_tags_.begin(); it != open_tags_.end(); ++it)
        {
          path_ += *it + ":";
        }
      }
      else if (element == "ITEMLIST")
      {
        std::vector<std::string> parts;

        if (list_.type == "string")
        {
          param_.setValue(list_.name, list_.stringlist, list_.description, list_.tags);
          if (list_.restrictions_index != -1)
          {
            StringUtils::split(list_.restrictions, ', ', parts);
            param_.setValidStrings(list_.name, ListUtils::create<std::string>(parts));
          }
        }
        else if (list_.type == "int")
        {
          param_.setValue(list_.name, list_.intlist, list_.description, list_.tags);
          if (list_.restrictions_index != -1)
          {
            StringUtils::split(list_.restrictions, ':', parts);
            if (parts.size() != 2)
            {
              StringUtils::split(list_.restrictions, '-', parts); //for downward compatibility
            }
            if (parts.size() == 2)
            {
              if (!parts[0].empty())
              {
                param_.setMinInt(list_.name, StringUtils::toInt32(parts[0]));
              }
              if (!parts[1].empty())
              {
                param_.setMaxInt(list_.name, StringUtils::toInt32(parts[1]));
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
            StringUtils::split(list_.restrictions, ':', parts);
            if (parts.size() != 2)
            {
              StringUtils::split(list_.restrictions, '-', parts); //for downward compatibility
            }
            if (parts.size() == 2)
            {
              if (!parts[0].empty())
              {
                param_.setMinFloat(list_.name, StringUtils::toDouble(parts[0]));
              }
              if (!parts[1].empty())
              {
                param_.setMaxFloat(list_.name, StringUtils::toDouble(parts[1]));
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
          warning(LOAD,StringUtils::toStr("Ignoring list entry '") + list_.name + "' because of unknown type '" + list_.type + "'");
        }
        list_.stringlist.clear();
        list_.intlist.clear();
        list_.doublelist.clear();
      }
    }

} // namespace OpenMS // namespace Internal
