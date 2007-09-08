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
      for (Int i=0; i<tag_num_; i++) in_tag_[i] = false;
      for (Int i=0; i<map_num_; i++) maps_[i] = Map();
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
      const XMLCh* xml_name = sm_.convert("name");
      const XMLCh* xml_value = sm_.convert("value");

      int tag = useMap_(TAGMAP,sm_.convert(qname),false,"opening tag");
      in_tag_[tag] = true;

      switch(tag)
      {
      case CELL:        cell_           = new GridCell(); break;
      case FPOSITION:   current_fcoord_ = asUInt_(sm_.convert(attributes.getValue(sm_.convert("dim")))); break;
      case SPOSITION:   current_scoord_ = asUInt_(sm_.convert(attributes.getValue(sm_.convert("dim")))); break;
      case PARAM:
        if (!(attributes.getIndex(xml_name)==-1) && !(attributes.getIndex(xml_value)==-1) )
        {
          param_->setValue(sm_.convert(attributes.getValue(xml_name)),sm_.convert(attributes.getValue(xml_value)));
        }
        break;
      case MAPPING:
        if (!(attributes.getIndex(xml_name)==-1))
        {
          String name = sm_.convert(attributes.getValue(xml_name));
          std::map<String,BaseMapping* >::const_iterator cit = mapping_instances_.find(name);
          if (cit == mapping_instances_.end())
          {
            const xercesc::Locator* loc = 0;
            setDocumentLocator(loc);
            String message = String("Error! This mapping type has not been registred with the XML Handler: ")+name;
            error(xercesc::SAXParseException(sm_.convert(message.c_str()), *loc));
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
      for (Int i=0; i<tag_num_; i++)
      {
        if (in_tag_[i])
        {
          GridCell::PositionType tmp;
          switch(i)
          {
          case FPOSITION:
            tmp = cell_->min();
            tmp[current_fcoord_] = asDouble_(sm_.convert(chars));
            cell_->setMin(tmp);
            break;
          case SPOSITION:
            tmp = cell_->max();
            tmp[current_scoord_] = asDouble_(sm_.convert(chars));
            cell_->setMax(tmp);
            break;
          }
        }
      }
    }

    // Docu in base class
    void GridHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      int tag = useMap_(TAGMAP,sm_.convert(qname),false,"closing tag");
      in_tag_[tag] = false;
      switch(tag)
      {
      case CELL:
        grid_->push_back(*cell_);
        delete cell_;
        break;
      case MAPPING:
        mapping_->setParameters(*param_);
        cell_->getMappings().push_back(mapping_);
        delete param_;
        registerMappings_();
        break;
      }
    }

    // Print the contents to a stream
    void GridHandler::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><!-- -*- mode: nxml; tab-width: 2 -*- -->" << std::endl;
      os << "<celllist>" << std::endl;

      // write features with their attributes
      for (UInt s=0; s<cgrid_->size(); s++)
      {
        const GridCell& cell = (*cgrid_)[s];

        os << "<cell nr=\"" << s << "\">" << std::endl;
        os << "\t<first>" << std::endl;
        DPosition<2> pos = cell.min();
        UInt dpos_size = pos.size();

        for (UInt i=0; i<dpos_size;i++)
        {
          os << "\t\t<fposition dim=\"" << i << "\">" << pos[i] << "</fposition>" <<  std::endl;
        }
        os << "\t</first>" << std::endl;

        os << "\t<second>" << std::endl;
        pos = cell.max();
        dpos_size = pos.size();

        for (UInt i=0; i<dpos_size;i++)
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
          Param map_param = (*citer)->getParameters();
          Param::ParamIterator piter = map_param.begin();
          while (piter != map_param.end())
          {
            os << "\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\">";
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






