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


#ifndef OPENMS_FORMAT_HANDLERS_FEATUREXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_FEATUREXMLHANDLER_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/TransService.hpp>

// STL includes
#include <iostream>


namespace OpenMS
{
	namespace Internal
	{

	/** 
		
		@brief XML Handler for a FeatureMap.
	 
		This class can be used to save the content of a
		FeatureMap into an XML file. The meta information
		(encapsulated by class ExperimentalSettings) is
		stored according to the mzData format. The features
		and their members are stored in a proprietary format
		(see funtion writeTo(stream& os) for details). 
	*/
  class FeatureXMLHandler
		: public SchemaHandler
  {
    public:
		typedef Feature::ConvexHullVector ConvexHullVector;
		typedef ConvexHullVector::value_type ConvexHullType;
						
      /**@name Constructors and destructor */
      //@{
      ///
      FeatureXMLHandler(FeatureMap<Feature>& map, const String& filename) 
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
			 	map_(&map), 
			 	cmap_(0),	
			 	feature_(), 
			 	exp_sett_()
  		{
				fillMaps_(Schemes::FeatureMap[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}
      
      ///
      FeatureXMLHandler(const FeatureMap<Feature>& map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				map_(0), 
				cmap_(&map),	
				feature_(), 
				exp_sett_()
  		{
				fillMaps_(Schemes::FeatureMap[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      ///
      virtual ~FeatureXMLHandler() 
      {
      }
      //@}

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, unsigned int length);

			///Writes the contents to a stream
			void writeTo(std::ostream& os);
			
			void setOptions(const PeakFileOptions& options)
			{ 
				options_ = options; 
			}

    protected:
		// Feature map pointer for reading
		FeatureMap<Feature>* map_;
		// Feature map pointer for writing
		const FeatureMap<Feature>* cmap_;

		/** @brief indices for tags used by FeatureXMLFile

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags { TAGNULL, FEATURELIST, FEATURE, POSITION, FEATINTENSITY, QUALITY, ACQUISITION,
								OVERALLQUALITY, CHARGE, FEATMODEL, PARAM, CONVEXHULL,
								HULLPOINT, HPOSITION, META, DESCRIPTION, FEATUREMAP, TAG_NUM};

		/** @brief indices for attributes used by FeatureXMLFile

			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Attributes { ATTNULL, DIM, NAME, VALUE, ATT_NUM};

		/** @brief indices for enum2str-maps used by FeatureXMLFile

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	TAGMAP, ATTMAP, MAP_NUM };

		PeakFileOptions options_;
		
		/**@name temporary datastructures to hold parsed data */
    //@{
		Feature* feature_;
		ModelDescription<2>* model_desc_;
		Param* param_;
		ConvexHullType* current_chull_;
		DPosition<2>* hull_position_;

		/// stream to collect experimental settings
		std::stringstream exp_sett_;
    //@}

 		// both quality and position might consist of several dimensions
 		// here we store the dimension that is currently parsed.
 		UInt current_pcoord_;				// current coordinate of the feature position
 		UInt current_qcoord_;				// coordinate of the feature quality
 		UInt current_hcoord_;				// coordinate of the current point in the hull

	};
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_FeatureXMLHandler_H
