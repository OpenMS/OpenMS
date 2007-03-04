// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/GridHandler.h>


namespace OpenMS
{
  namespace Internal
  {
    GridHandler::GridHandler(Grid& grid, const String& filename)
        : XMLHandler(filename),
        grid_(&grid),
        cgrid_(0),
        cell_(),
        mapping_(),
        param_()
    {
      for (SignedInt i=0; i<TAG_NUM; i++) in_tag_[i] = false;
      for (SignedInt i=0; i<MAP_NUM; i++) maps[i] = Map();
      setConstants_();
      fillMaps_();
      registerMappings_();
    }

    GridHandler::GridHandler(const Grid& grid, const String& filename)
        : XMLHandler(filename),
        grid_(0),
        cgrid_(&grid),
        cell_(),
        mapping_(),
        param_()
    {
      setConstants_();
      fillMaps_();
      registerMappings_();
    }

    void GridHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      const XMLCh* xml_name = xercesc::XMLString::transcode("name");
      const XMLCh* xml_value = xercesc::XMLString::transcode("value");

      int tag = useMap_(TAGMAP,xercesc::XMLString::transcode(qname),false,"opening tag");
      in_tag_[tag] = true;

      switch(tag)
      {
      case CELL:        cell_           = new GridCell(); break;
      case FPOSITION:   current_fcoord_ = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("dim")))); break;
      case SPOSITION:   current_scoord_ = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("dim")))); break;
      case PARAM:
        if (!(attributes.getIndex(xml_name)==-1) && !(attributes.getIndex(xml_value)==-1) )
        {
          param_->setValue(xercesc::XMLString::transcode(attributes.getValue(xml_name)),xercesc::XMLString::transcode(attributes.getValue(xml_value)));
        }
        break;
      case MAPPING:
        if (!(attributes.getIndex(xml_name)==-1))
        {
          String name = xercesc::XMLString::transcode(attributes.getValue(xml_name));
          std::map<String,BaseMapping* >::const_iterator cit = mapping_instances.find(name);
          if (cit == mapping_instances.end())
          {
            const xercesc::Locator* loc = 0;
            setDocumentLocator(loc);
            String message = String("Error! This mapping type has not been registred with the XML Handler: ")+name;
            error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc));
          }
          else
          {
            param_   = new Param();
            mapping_ = cit->second;
          }
        } // end if (!attributes..)
        break;
      }
    }

    // Docu in base class
    void GridHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
    {
      for (SignedInt i=0; i<TAG_NUM; i++)
      {
        if (in_tag_[i])
        {
          GridCell::PositionType tmp;
          switch(i)
          {
          case FPOSITION:
            tmp = cell_->min();
            tmp[current_fcoord_] = asDouble_(xercesc::XMLString::transcode(chars));
            cell_->setMin(tmp);
            break;
          case SPOSITION:
            tmp = cell_->max();
            tmp[current_scoord_] = asDouble_(xercesc::XMLString::transcode(chars));
            cell_->setMax(tmp);
            break;
          }
        }
      }
    }

    // Docu in base class
    void GridHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      int tag = useMap_(TAGMAP,xercesc::XMLString::transcode(qname),false,"closing tag");
      in_tag_[tag] = false;
      switch(tag)
      {
      case CELL:
        grid_->push_back(*cell_);
        delete cell_;
        break;
      case MAPPING:
        mapping_->setParam(*param_);
        cell_->getMappings().push_back(mapping_);
        delete param_;
        registerMappings_();
        break;
      }
    }

    /// Print the contents to a stream
    void GridHandler::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><!-- -*- mode: nxml; tab-width: 2 -*- -->" << std::endl;
      os << "<celllist>" << std::endl;

      // write features with their attributes
      for (UnsignedInt s=0; s<cgrid_->size(); s++)
      {
        const GridCell& cell = (*cgrid_)[s];

        os << "<cell nr=\"" << s << "\">" << std::endl;
        os << "\t<first>" << std::endl;
        DPosition<2> pos = cell.min();
        UnsignedInt dpos_size = pos.size();

        for (UnsignedInt i=0; i<dpos_size;i++)
        {
          os << "\t\t<fposition dim=\"" << i << "\">" << pos[i] << "</fposition>" <<  std::endl;
        }
        os << "\t</first>" << std::endl;

        os << "\t<second>" << std::endl;
        pos = cell.max();
        dpos_size = pos.size();

        for (UnsignedInt i=0; i<dpos_size;i++)
        {
          os << "\t\t<sposition dim=\"" << i << "\">" << pos[i] << "</sposition>" <<  std::endl;
        }
        os << "\t</second>" << std::endl;


        os << "\t<mappinglist>" << std::endl;
        MappingVector mappings = cell.getMappings();

        MappingVector::const_iterator citer = mappings.begin();

        while (citer != mappings.end() )
        {
          os << "\t\t<mapping name=\"" << (*citer)->getName() << "\">" << std::endl;
          Param map_param = (*citer)->getParam();
          Param::ConstIterator piter = map_param.begin();
          while (piter != map_param.end())
          {
            os << "\t\t\t<param name=\"" << piter->first << "\" value=\"" << piter->second << "\">";
            os << "</param>" << std::endl;
            piter++;
          }
          os << "\t\t</mapping>" << std::endl;
          citer++;
        }

        os << "\t</mappinglist>" << std::endl;
        os << "</cell>" << std::endl;

      } // end for ( features )

      os << "</celllist>" << std::endl;
    }
  } // namespace Internal
} // namespace OpenMS






