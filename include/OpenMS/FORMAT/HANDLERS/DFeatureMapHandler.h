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


#ifndef OPENMS_FORMAT_HANDLERS_DFEATUREMAPHANDLER_H
#define OPENMS_FORMAT_HANDLERS_DFEATUREMAPHANDLER_H

#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/FORMAT/Param.h>
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
		
		@brief XML Handler for a DFeatureMap.
	 
		This class can be used to save the content of a
		DFeatureMap into an XML file. The meta information
		(encapsulated by class ExperimentalSettings) is
		stored according to the mzData format. The features
		and their members are stored in a proprietary format
		inspired by mzData (see funtion writeTo(stream& os)
		for details. 
	*/
  template <Size D, typename TraitsT = KernelTraits, typename FeatureT = DFeature<D> >
  class DFeatureMapHandler
		: public SchemaHandler
  {
    public:
    	/**
				@name Type definitions
			*/
			//@{
			typedef TraitsT TraitsType;
			typedef FeatureT FeatureType;
			typedef std::vector<FeatureType> ContainerType;
			typedef typename ContainerType::iterator Iterator;
			typedef typename ContainerType::const_iterator ConstIterator;
			typedef typename ContainerType::reverse_iterator ReverseIterator;
			typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
			typedef FeatureType& Reference;
			typedef const FeatureType& ConstReference;
			typedef typename FeatureType::ConvexHullVector ConvexHullVector;
			typedef typename FeatureType::ConvexHullType ConvexHullType;

			// STL compatibility
			typedef FeatureType value_type;
			typedef FeatureType* pointer;
			typedef const FeatureType* const_pointer;
			typedef Reference reference;
			typedef ConstReference const_reference;
			typedef typename ContainerType::size_type size_type;
			typedef typename ContainerType::difference_type difference_type;
			typedef Iterator iterator;
			typedef ConstIterator const_iterator;
			typedef ReverseIterator reverse_iterator;
			typedef ConstReverseIterator const_reverse_iterator;
			//@}
						
      /**@name Constructors and destructor */
      //@{
      ///
      DFeatureMapHandler(DFeatureMap<D,FeatureType>& map, const String& filename) 
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
			 	map_(&map), 
			 	cmap_(0),	
			 	feature_(), 
			 	exp_sett_()
  		{
				fillMaps_(Schemes::DFeatureMap[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}
      
      ///
      DFeatureMapHandler(const DFeatureMap<D,FeatureType>& map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				map_(0), 
				cmap_(&map),	
				feature_(), 
				exp_sett_()
  		{
				fillMaps_(Schemes::DFeatureMap[schema_]);	// fill maps with current schema
				setMaps_(TAGMAP, ATTMAP);
			}

      ///
      virtual ~DFeatureMapHandler() { }
      //@}

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);

			///Writes the contents to a stream
			void writeTo(std::ostream& os);
			
			void setOptions(const PeakFileOptions& options) { options_ = options; }

    protected:
		// Feature map pointer for reading
		DFeatureMap<D, FeatureType>* map_;
		// Feature map pointer for writing
		const DFeatureMap<D, FeatureType>* cmap_;

		/** @brief indices for tags used by DFeatureMapFile

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags { TAGNULL, FEATURELIST, FEATURE, POSITION, FEATINTENSITY, QUALITY, ACQUISITION,
								OVERALLQUALITY, CHARGE, FEATMODEL, PARAM, CONVEXHULL,
								HULLPOINT, HPOSITION, META, DESCRIPTION, FEATUREMAP, TAG_NUM};

		/** @brief indices for attributes used by DFeatureMapFile

			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Attributes { ATTNULL, DIM, NAME, VALUE, ATT_NUM};

		/** @brief indices for enum2str-maps used by DFeatureMapFile

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	TAGMAP, ATTMAP, MAP_NUM };

		PeakFileOptions options_;
		
		/**@name temporary datastructures to hold parsed data */
    //@{
		FeatureType* feature_;
		ModelDescription<D>* model_desc_;
		Param* param_;
		ConvexHullType* current_chull_;
		DPosition<D>* hull_position_;

		/// stream to collect experimental settings
		std::stringstream exp_sett_;
    //@}

 		// both quality and position might consist of several dimensions
 		// here we store the dimension that is currently parsed.
 		Position current_pcoord_;				// current coordinate of the feature position
 		Position current_qcoord_;				// coordinate of the feature quality
 		Position current_hcoord_;				// coordinate of the current point in the hull

	};



	//--------------------------------------------------------------------------------

  template <Size D, typename TraitsT, typename FeatureT>
  void DFeatureMapHandler<D,TraitsT,FeatureT>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << "</" << xercesc::XMLString::transcode(qname) << ">\n";
			if (String(xercesc::XMLString::transcode(qname)) != enum2str_(TAGMAP,DESCRIPTION))
			{
				return;
			}
		}

		int tag = leaveTag(qname);

		// Do something depending on the tag
		switch(tag) {
			case DESCRIPTION:
				// delegate control to ExperimentalSettings handler
				{
					// initialize parser
					xercesc::XMLPlatformUtils::Initialize();
					xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
					parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
					parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
					
					MzDataExpSettHandler handler( *((ExperimentalSettings*)map_),file_);
					handler.resetErrors();
					parser->setContentHandler(&handler);
					parser->setErrorHandler(&handler);
					
					String tmp(exp_sett_.str().c_str());
					
					xercesc::MemBufInputSource source((const XMLByte*)(tmp.c_str()), tmp.size(), "dummy");
	      	parser->parse(source);
	      	delete(parser);
				}
				break;
			case FEATURE:
				map_->push_back(*feature_);
				delete feature_;
				break;
			case FEATMODEL:
				model_desc_->setParam(*param_);
				feature_->setModelDescription(*model_desc_);
				delete param_;
				delete model_desc_;
				break;
			case HULLPOINT:
				current_chull_->addPoint(*hull_position_);
				delete hull_position_;
				break;
			case CONVEXHULL:
				feature_->getConvexHulls().push_back(*current_chull_);
				delete current_chull_;
				break;
		}
  }

  template <Size D, typename TraitsT, typename FeatureT>
  void DFeatureMapHandler<D,TraitsT,FeatureT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << '<' << xercesc::XMLString::transcode(qname);
			Size n=attributes.getLength();
			for (Size i=0; i<n; ++i)
			{
				exp_sett_ << ' ' << xercesc::XMLString::transcode(attributes.getQName(i)) << "=\""	<< xercesc::XMLString::transcode(attributes.getValue(i)) << '\"';
			}
			exp_sett_ << '>';
			return;
		}

		int tag = enterTag(qname, attributes);

		String tmp_str;
		// Do something depending on the tag
		switch(tag) {
			case DESCRIPTION: 
				exp_sett_ << '<' << xercesc::XMLString::transcode(qname) << '>'; 
				break;
   		case FEATURE: 	 
   			feature_        = new DFeature<D>();
   			break;
			case QUALITY:
				tmp_str = getAttributeAsString_(DIM);
				current_qcoord_ = asUnsignedInt_(tmp_str); 
				break;
			case POSITION:
				tmp_str = getAttributeAsString_(DIM);
				current_pcoord_ = asUnsignedInt_(tmp_str); 
				break;
			case CONVEXHULL: 
				current_chull_  = new ConvexHullType(); 
				break;
			case HULLPOINT:  
				hull_position_  = new DPosition<D>(); 
				break;
			case HPOSITION:  
				tmp_str = getAttributeAsString_(DIM);
				current_hcoord_ = asUnsignedInt_(tmp_str); 
				break;
			case FEATMODEL:
				model_desc_ = new ModelDescription<D>();
		  	param_ = new Param();
				tmp_str = getAttributeAsString_(NAME);
		  	if (tmp_str != "")
		  	{
		  		model_desc_->setName(tmp_str);
		  	}
		  	break;
		  case PARAM:
		  {
		  	String name = getAttributeAsString_(NAME);
				String value = getAttributeAsString_(VALUE);
		  	if (name != "" && value != "")
		  		param_->setValue(name, value);
		  	break;
		  }
		}
	}

	template <Size D, typename TraitsT, typename FeatureT>
  void DFeatureMapHandler<D,TraitsT,FeatureT>::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << xercesc::XMLString::transcode(chars);
			return;
		}

		// find the tag that the parser is in right now
 		for (Size i=0; i<is_parser_in_tag_.size(); i++)
 		{
			if (is_parser_in_tag_[i])
			{
				switch(i) 
				{
					case FEATINTENSITY: 
						feature_->setIntensity(asDouble_(xercesc::XMLString::transcode(chars))); 
						break;
					case POSITION:
						feature_->getPosition()[current_pcoord_] = asDouble_(xercesc::XMLString::transcode(chars));
						break;
					case QUALITY:       
						feature_->getQuality(current_qcoord_) = asDouble_(xercesc::XMLString::transcode(chars));
							break;
					case OVERALLQUALITY:  
						feature_->getOverallQuality() = asDouble_(xercesc::XMLString::transcode(chars)); break;
					case CHARGE:          
						feature_->setCharge(asSignedInt_(xercesc::XMLString::transcode(chars)));
						break;
					case HPOSITION:       
						(*hull_position_)[current_hcoord_] = asDouble_(xercesc::XMLString::transcode(chars)); 
						break;
					case META:						
						feature_->setMetaValue(3,String(xercesc::XMLString::transcode(chars))); 
						break;
				}
			}
		}
  }


	template <Size D, typename TraitsT, typename FeatureT>
 	void DFeatureMapHandler<D,TraitsT,FeatureT>::writeTo(std::ostream& os)
	{
		UniqueIdGenerator id_generator = UniqueIdGenerator::instance();

		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
		   << "<featureMap>\n";

		// delegate control to ExperimentalSettings handler
		Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cmap_),"");
		handler.writeTo(os);

		os << "\t<featureList count=\"" << cmap_->size() << "\">\n";

		// write features with their corresponding attributes
		for (UnsignedInt s=0; s<cmap_->size(); s++)
		{
			const DFeature<D>& dfeat = (*cmap_)[s];

			os << "\t\t<feature id=\"" << id_generator.getUID() << "\">" << std::endl;

			DPosition<D> pos = dfeat.getPosition();
			UnsignedInt dpos_size = pos.size();
			for (UnsignedInt i=0; i<dpos_size;i++)
				os <<	"\t\t\t<position dim=\"" << i << "\">" << pos[i] << "</position>" << 	std::endl;

			os << "\t\t\t<intensity>" << dfeat.getIntensity() << "</intensity>" << std::endl;

			for (UnsignedInt i=0; i<dpos_size;i++)
			os << "\t\t\t<quality dim=\"" << i << "\">" << dfeat.getQuality(i) << "</quality>" << std:: endl;

			if(dfeat.getMetaValue(3)!=DataValue::EMPTY)
				os << "\t\t\t<meta>" << dfeat.getMetaValue(3) << "</meta>" << std:: endl;

			os << "\t\t\t<overallquality>" << dfeat.getOverallQuality() << "</overallquality>" << std:: endl;
			os << "\t\t\t<charge>" << dfeat.getCharge() << "</charge>" << std:: endl;

			// write model description
			ModelDescription<D> desc = dfeat.getModelDescription();
			os << "\t\t\t<model name=\"" << desc.getName() << "\">" << std:: endl;
			Param modelp = desc.getParam();
			Param::ConstIterator piter = modelp.begin();
			while (piter != modelp.end())
			{
				os << "\t\t\t\t<param name=\"" << piter->first << "\" value=\"" << piter->second << "\">";
				os << "</param>" << std::endl;
				piter++;
			}
			os << "\t\t\t</model>" << std::endl;

			// write convex hull
			DFeature<2>::ConvexHullVector hulls = dfeat.getConvexHulls();
			DFeature<2>::ConvexHullVector::iterator citer = hulls.begin();

			UnsignedInt hulls_count = hulls.size();

			for (UnsignedInt i=0;i<hulls_count; i++)
			{
				os << "\t\t\t<convexhull nr=\"" << i << "\">" << std:: endl;

				ConvexHullType current_hull = hulls[i];
				UnsignedInt hull_size       = current_hull.getPoints().size();

				for (UnsignedInt j=0;j<hull_size;j++)
				{
					os << "\t\t\t\t<hullpoint>" << std::endl;

					DPosition<D> pos = current_hull.getPoints()[j];
					UnsignedInt pos_size = pos.size();
					for (UnsignedInt k=0; k<pos_size; k++)
					{
						os << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << std::endl;
					}

					os << "\t\t\t\t</hullpoint>" << std::endl;
				} // end for (..hull_size..)

				os << "\t\t\t</convexhull>" << std::endl;
			} // end  for ( ... hull_count..)

			os << "\t\t</feature>\n";

		} // end for ( features )

		os << "\t</featureList>\n</featureMap>\n";
		os <<
			"<!-- Local Variables: -->\n"
			"<!-- mode: nxml -->\n"
			"<!-- tab-width: 2 -->\n"
			"<!-- End: -->\n";
	}


	} // namespace Internal
} // namespace OpenMS

#endif


// Please leave the page-feed character (Ctrl-L) in the following line,
// otherwise emacs will not select c++ mode when opening this file.

// Thanks!

// EOF
