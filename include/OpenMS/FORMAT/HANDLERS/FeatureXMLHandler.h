// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
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
		 
			This class can be used to save the content of a FeatureMap into an XML file. The meta information
			(encapsulated by class ExperimentalSettings) is stored according to the mzData format. The features
			and their members are stored in a proprietary format (see funtion writeTo(stream& os) for details). 
		*/
	  class FeatureXMLHandler
			: public XMLHandler
	  {
	    public:
	      /**@name Constructors and destructor */
	      //@{
	      /// Constructor for reading 
	      FeatureXMLHandler(FeatureMap<Feature>& map, const String& filename, const String& version) 
	      : XMLHandler(filename,version),
				 	map_(&map), 
				 	cmap_(0),	
				 	exp_sett_(),
				 	in_description_(false)
	  		{
				}
	      
	      ///Constructor for writing
	      FeatureXMLHandler(const FeatureMap<Feature>& map, const String& filename, const String& version)
	      : XMLHandler(filename,version),
					map_(0), 
					cmap_(&map),	
					exp_sett_(),
				 	in_description_(false)
	  		{
				}
	
	      /// Destructor
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
				
				///Sets the options
				void setOptions(const PeakFileOptions& options);
	
	    protected:
				/// Feature map pointer for reading
				FeatureMap<Feature>* map_;
				/// Feature map pointer for writing
				const FeatureMap<Feature>* cmap_;
				/// Options that can be set				
				PeakFileOptions options_;
				
				/**@name temporary datastructures to hold parsed data */
		    //@{
				Feature feature_;
				ModelDescription<2>* model_desc_;
				Param param_;
				ConvexHull2D current_chull_;
				DPosition<2> hull_position_;	
				/// stream to collect experimental settings
				std::stringstream exp_sett_;
		    //@}
				
				/// current dimension of the feature position, quality, or convex hull point
		 		UInt dim_;			
				
				//flag that indicates that the parser in in the description secion
				bool in_description_;
		};
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_FeatureXMLHandler_H
