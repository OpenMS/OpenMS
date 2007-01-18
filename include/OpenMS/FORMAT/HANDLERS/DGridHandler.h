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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_DGRIDHANDLER_H
#define OPENMS_FORMAT_HANDLERS_DGRIDHANDLER_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

// all implementations of class DBaseMapping must be
// included here 
#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>

// STL includes
#include <iostream>
#include <valarray>
#include <string>

#include <xercesc/sax2/Attributes.hpp>

namespace OpenMS
{
	namespace Internal
	{

	/** @brief XML Handler for a vector of grid cells including their transformations.
	  
	  	This is a simplified version of class DFeatureMapHandler. We explicitly allow
	  	several tagtypes even if just one type is used in this implementation (for
	  	details see class member further below). Therefore this class can be extended
	  	in the future in order to save meta information with the grid such as information
	  	about the experiment etc.
	  	
	  	@note A grid cell can have different transformations for each dimension.
	  	If you want this XML handler class to support other transformations than the
	  	linear one, you must register this class with the handler. For details, have a look
	  	at registerMappings_() .
	 */
  template <Size D>
  class DGridHandler
		: public XMLHandler
  {
    public:
    	/**	
				@name Type definitions
			*/
			//@{
			typedef typename DGridCell<D>::MappingVector MappingVector;
			//@}
    							
      /**@name Constructors and destructor */
      //@{
      ///
      DGridHandler(DGrid<D>& grid, const String& filename) 
      : XMLHandler(filename),
      	grid_(&grid), 
      	cgrid_(0), 
				cell_(), 
				mapping_(), 
				param_()
  		{
				for (Index i=0; i<TAG_NUM; i++)	in_tag_[i] = false;
				for (Index i=0; i<MAP_NUM; i++)	maps[i] = Map();
				setConstants_();
				fillMaps_();
				registerMappings_();
			}
      
      ///
      DGridHandler(const DGrid<D>& grid, const String& filename)
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
      ///
      virtual ~DGridHandler()  
      {
      }     
      //@}

			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  		{
  			
  			const XMLCh* xml_name = xercesc::XMLString::transcode("name");
  			const XMLCh* xml_value = xercesc::XMLString::transcode("value");
  			
  			int tag = useMap_(TAGMAP,xercesc::XMLString::transcode(qname),false,"opening tag");
				in_tag_[tag] = true;
				
				switch(tag) 
				{
					case CELL: 				cell_           = new DGridCell<D>(); break;
					case FPOSITION:		current_fcoord_ = asUnsignedInt_(xercesc::XMLString::transcode(attributes.getValue(xercesc::XMLString::transcode("dim")))); break;
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
		  				typename std::map<String,DBaseMapping<1>* >::const_iterator cit = mapping_instances.find(name);
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
      virtual void characters(const XMLCh* const chars, const unsigned int /*length*/)
      {
      	for (Index i=0; i<TAG_NUM; i++) 
      	{
						if (in_tag_[i])
						{
							typename DGridCell<D>::PositionType tmp;
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
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
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
			void writeTo(std::ostream& os)
			{
					
				os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><!-- -*- mode: nxml; tab-width: 2 -*- -->" << std::endl;
				os << "<celllist>" << std::endl;	
						 
				// write features with their attributes				
				for (UnsignedInt s=0; s<cgrid_->size(); s++)
				{
					const DGridCell<D>& cell = (*cgrid_)[s];
					
					os << "<cell nr=\"" << s << "\">" << std::endl;
					os << "\t<first>" << std::endl;
					DPosition<D> pos = cell.min();
					UnsignedInt dpos_size = pos.size();

					for (UnsignedInt i=0; i<dpos_size;i++)
					{
						os <<	"\t\t<fposition dim=\"" << i << "\">" << pos[i] << "</fposition>" << 	std::endl; 	
					}
					os << "\t</first>" << std::endl;
					
					os << "\t<second>" << std::endl;
					pos = cell.max();
					dpos_size = pos.size();
			
					for (UnsignedInt i=0; i<dpos_size;i++)
					{
						os <<	"\t\t<sposition dim=\"" << i << "\">" << pos[i] << "</sposition>" << 	std::endl; 	
					}
					os << "\t</second>" << std::endl;
					
					
					os << "\t<mappinglist>" << std::endl;
					MappingVector mappings = cell.getMappings();

					typename MappingVector::const_iterator citer = mappings.begin();
					
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
		
		protected:
		
		std::vector<String> tagsVector_;	
				
		/// Maps to assoziate Strings with enumeration values
		enum MapType {	TAGMAP };
		typedef std::map<std::string,int> Map;
		static const int MAP_NUM = 1;
		Map maps[MAP_NUM];
		
		/// Vector of grid cell to be read
		DGrid<D>* grid_;
		/// Vector of pairs to be written
		const DGrid<D>* cgrid_;
		
		/// The tags we expect to encounter
		enum Tags { CELLLIST, CELL, FIRSTPOSITION, SECONDPOSITION, 
			          FPOSITION, SPOSITION, MAPPINGLIST, MAPPING, PARAM };
						
		static const int TAG_NUM = 9;
					
		/// Indicates which tag is currently parsed
		bool in_tag_[TAG_NUM];

		// temporary datastructures to hold parsed data
		DGridCell<D>* cell_;
		DBaseMapping<1>* mapping_;
		Param* param_;

		Position current_fcoord_;				
		Position current_scoord_;
		
		std::map<String,DBaseMapping<1>* > mapping_instances;
					
		void fillMaps_() 
		{
			fillMap_(maps[TAGMAP], tagsVector_);
		}
		
		/// mapping types must be registred with the handler class
		void registerMappings_()
		{
			// insert new mappings (transformations) here.
			mapping_instances["DLinearMapping"] = new DLinearMapping<1>();
		}
		
		/// @brief Find value in the given map
		/// if not found: fatal error or warning message
		inline int useMap_(MapType type, String value, bool fatal=true, const char* message="")
		{
			Map::const_iterator it =  maps[type].find(value);
			if (it == maps[type].end())
			{
				if (fatal)
				{
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("Error in enumerated value \"") + value + "\"";
					error(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
				}
				else if (message != "")
				{
					const xercesc::Locator* loc = 0;
					setDocumentLocator(loc);
					String message = String("Unhandled ") + message + "\"" + value + "\"";
					warning(xercesc::SAXParseException(xercesc::XMLString::transcode(message.c_str()), *loc ));
				}
			}	
			else 
			{
				return it->second;
			}
			return 0;
		}

		///  @brief Create map from the given set of strings
		inline void fillMap_(Map& dict, std::vector<String> array)
		{
			for (UnsignedInt i=0; i<array.size(); i++)
			{
				 dict[ std::string(array[i]) ] = i; 
			}
		}
		
		/** @brief Set constants of XML handler */
		inline void setConstants_()
		{ 					 	
			char* tags[] 
			= {"celllist", "cell", "first", "second", "fposition", 
				 "sposition", "mappinglist", "mapping", "param"};
			fillVector_(tagsVector_,tags,TAG_NUM);
		}
		
		inline void fillVector_(std::vector<String>& vec, char* contents[], int nr)
		{
			for (int i=0; i<nr;i++)	vec.push_back(contents[i]);					
		}
				 		
  }; // end of class DGridHandler

	} // namespace Internal
} // namespace OpenMS

#endif
