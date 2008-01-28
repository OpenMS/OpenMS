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
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Grid.h>

namespace OpenMS
{
  namespace Internal
  {
    GridHandler::GridHandler(Grid& grid, const String& filename)
    	: XMLHandler(filename),
        grid_(&grid),
        cgrid_(0),
        mapping_(0),
        param_()
    {
    }

    GridHandler::GridHandler(const Grid& grid, const String& filename)
    	: XMLHandler(filename),
        grid_(0),
        cgrid_(&grid),
        mapping_(0),
        param_()
    {
    }

    GridHandler::~GridHandler()
    {
    }
			
    void GridHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      const XMLCh* s_name = xercesc::XMLString::transcode("name");
      const XMLCh* s_value = xercesc::XMLString::transcode("value");
      const XMLCh* s_dim = xercesc::XMLString::transcode("dim");
      
      String tag = sm_.convert(qname);
			open_tags_.push_back(tag);

      if (tag=="cell")
      {
      	grid_->resize(grid_->size()+1);
      }
      else if (tag=="fposition")
      {
      	dim_ = attributeAsInt_(attributes, s_dim);
      }
      else if (tag=="sposition")
      {
      	dim_ = attributeAsInt_(attributes, s_dim);
      }
      else if (tag=="param")
      {
      	param_.setValue(attributeAsString_(attributes, s_name),attributeAsString_(attributes, s_value));
      }
      else if (tag=="mapping")
      {
      	String name = attributeAsString_(attributes, s_name);
				if (name=="LinearMapping")
				{
          param_ = Param();
          mapping_ = new LinearMapping();
				}
        else
        {
					error(String("Error! This mapping type has not been registered with the XML Handler: ")+name, 0, 0);
        }
      }
    }

    // Docu in base class
    void GridHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
    {
			if (open_tags_.back()=="fposition")
			{
        GridCell::PositionType tmp = grid_->back().min();
        tmp[dim_] = asDouble_(sm_.convert(chars));
        grid_->back().setMin(tmp);
      }
      else if (open_tags_.back()=="sposition")
			{
  			GridCell::PositionType tmp = grid_->back().max();
        tmp[dim_] = asDouble_(sm_.convert(chars));
        grid_->back().setMax(tmp);
      }
    }

    // Docu in base class
    void GridHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
    	const XMLCh* s_mapping = xercesc::XMLString::transcode("mapping");
    	
			if (equal_(qname,s_mapping))
	    {
        mapping_->setParameters(param_);
        param_ = Param();
        grid_->back().getMappings().push_back(mapping_);
      }
    	
    	open_tags_.pop_back();
    }

    // Print the contents to a stream
    void GridHandler::writeTo(std::ostream& os)
    {
      os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
      os << "<celllist>" << std::endl;

      // write features with their attributes
      for (UInt s=0; s<cgrid_->size(); s++)
      {
        const GridCell& cell = (*cgrid_)[s];

        os << "\t<cell nr=\"" << s << "\">" << std::endl;
        os << "\t\t<first>" << std::endl;
        DPosition<2> pos = cell.min();
        UInt dpos_size = pos.size();

        for (UInt i=0; i<dpos_size;i++)
        {
          os << "\t\t\t<fposition dim=\"" << i << "\">" << pos[i] << "</fposition>" <<  std::endl;
        }
        os << "\t\t</first>" << std::endl;

        os << "\t\t<second>" << std::endl;
        pos = cell.max();
        dpos_size = pos.size();

        for (UInt i=0; i<dpos_size;i++)
        {
          os << "\t\t\t<sposition dim=\"" << i << "\">" << pos[i] << "</sposition>" <<  std::endl;
        }
        os << "\t\t</second>" << std::endl;


        os << "\t\t<mappinglist>" << std::endl;
        GridCell::MappingVector mappings = cell.getMappings();

        GridCell::MappingVector::const_iterator citer = mappings.begin();

        while (citer != mappings.end() )
        {
          os << "\t\t\t<mapping name=\"" << (*citer)->getName() << "\">" << std::endl;
          Param map_param = (*citer)->getParameters();
          Param::ParamIterator piter = map_param.begin();
          while (piter != map_param.end())
          {
            os << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\">" << "</param>" << std::endl;
            piter++;
          }
          os << "\t\t\t</mapping>" << std::endl;
          citer++;
        }

        os << "\t\t</mappinglist>" << std::endl;
        os << "\t</cell>" << std::endl;

      } // end for ( features )

      os << "</celllist>" << std::endl;
    }
  } // namespace Internal
} // namespace OpenMS
