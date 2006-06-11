// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DFeatureMapHandler.h,v 1.17 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#ifndef OPENMS_FORMAT_HANDLERS_DFEATUREMAPHANDLER_H
#define OPENMS_FORMAT_HANDLERS_DFEATUREMAPHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/FORMAT/Param.h>

// STL includes
#include <iostream>
#include <valarray>
#include <string>

namespace OpenMS
{
	namespace Internal
	{

	/** @brief XML Handler for a DFeatureMap.
	 * 
	 *  This class can be used to save the content of a
	 *  DFeatureMap into an XML file. The meta information
	 *  (encapsulated by class ExperimentalSettings) is
	 *  stored according to the mzData format. The features
	 *  and their members are stored in a proprietary format
	 *  inspired by mzData (see funtion writeTo(stream& os)
	 *  for details.
	 * 
	 * */
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
      DFeatureMapHandler(DFeatureMap<D,TraitsType,FeatureType>& map) 
      : SchemaHandler(TAG_NUM,MAP_NUM),
			 	map_(&map), cmap_(0),	feature_(), exp_sett_(exp_sett_str_, IO_ReadWrite)
  		{
				file_ = __FILE__;
				fillMaps_(Schemes::DFeatureMap[schema_]);	// fill maps with current schema
			}
      
      ///
      DFeatureMapHandler(const DFeatureMap<D,TraitsType,FeatureType>& map)
      : SchemaHandler(TAG_NUM,MAP_NUM),
				map_(0), cmap_(&map),	feature_(), exp_sett_(exp_sett_str_, IO_ReadWrite)
  		{
				file_ = __FILE__;
				fillMaps_(Schemes::DFeatureMap[schema_]);	// fill maps with current schema
			}

      ///
      virtual ~DFeatureMapHandler() { }
      //@}
			
			/// Show warnings of the QT parser or not.
			inline void showWarnings(bool w) { use_warnings_ = w; }

      /// This function is called for each closing tag in the XML file.
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname );

      /// This function is called for each opening XML tag in the file.
      virtual bool startElement(const QString & uri, const QString & local_name,
																const QString & qname, const QXmlAttributes & attributes );

			/// This function is called by the parser for each chunk of
		  /// characters between two tags.
      virtual bool characters( const QString & chars );

			///Writes the contents to a stream
			void writeTo(std::ostream& os);

    protected:
		// Feature map pointer for reading
		DFeatureMap<D,TraitsType, FeatureType>* map_;
		// Feature map pointer for writing
		const DFeatureMap<D,TraitsType, FeatureType>* cmap_;

		/** @brief indices for tags used by DFeatureMapFile

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags { TAGNULL, FEATURELIST, FEATURE, POSITION, FEATINTENSITY, QUALITY, ACQUISITION,
								OVERALLQUALITY, CHARGE, FEATMODEL, PARAM, CONVEXHULL,
								HULLPOINT, HPOSITION, META, DESCRIPTION, FEATUREMAP, TAG_NUM};

		/** @brief indices for enum2str-maps used by DFeatureMapFile

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	TAGMAP, MAP_NUM };

		/**@name temporary datastructures to hold parsed data */
    //@{
		FeatureType* feature_;
		ModelDescription<D>* model_desc_;
		Param* param_;
		ConvexHullType* current_chull_;
		DPosition<D>* hull_position_;

		/// stream to collect experimental settings
		QTextStream exp_sett_;
		/// string with the xml containing the experimental settings
		QString exp_sett_str_;
    //@}

 		// both quality and position might consist of several dimensions
 		// here we store the dimension that is currently parsed.
 		Position current_pcoord_;				// current coordinate of the feature position
 		Position current_qcoord_;				// coordinate of the feature quality
 		Position current_hcoord_;				// coordinate of the current point in the hull

	};



	//--------------------------------------------------------------------------------

  template <Size D, typename TraitsT, typename FeatureT>
  bool DFeatureMapHandler<D,TraitsT,FeatureT>::endElement( const QString & /*uri*/,
	const QString & /*local_name*/,	 const QString & qname )
	{
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << "</" << qname << ">\n";
			if (qname != enum2str_(TAGMAP,DESCRIPTION))	return true;
		}

	  int tag = str2enum_(TAGMAP,qname,"closing tag");  // index of current tag
		is_parser_in_tag_[tag] = false;

		// Do something depending on the tag
		switch(tag) {
			case DESCRIPTION:
				// delegate control to ExperimentalSettings handler
				{
					QXmlSimpleReader parser;
					srand(static_cast<unsigned>(time(0)));
					parser.setFeature("http://xml.org/sax/features/namespaces",false);
					parser.setFeature("http://xml.org/sax/features/namespace-prefixes", false);

					MzDataExpSettHandler handler( *((ExperimentalSettings*)map_));
					parser.setContentHandler(&handler);
					parser.setErrorHandler(&handler);
					QXmlInputSource source;
					source.setData(exp_sett_str_);
					parser.parse(source);
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
				current_chull_->push_back(*hull_position_);
				delete hull_position_;
				break;
			case CONVEXHULL:
				feature_->getConvexHulls().push_back(*current_chull_);
				delete current_chull_;
				break;
		}
		return true;
  }

  template <Size D, typename TraitsT, typename FeatureT>
  bool DFeatureMapHandler<D,TraitsT,FeatureT>::startElement(
				const QString & /*uri*/, const QString & /*local_name*/,
				const QString & qname, const QXmlAttributes & attributes)
	{
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << '<' << qname;
			Size n=attributes.count();
			for (Size i=0; i<n; ++i)
				exp_sett_ << ' ' << attributes.qName(i) << "=\""	<< attributes.value(i) << '\"';
			exp_sett_ << '>';
			return true;
		}

		int tag = str2enum_(TAGMAP,qname,"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;

		// Do something depending on the tag
		switch(tag) {
			case DESCRIPTION: exp_sett_ << '<' << qname << '>'; break;
   		case FEATURE: 	 feature_        = new DFeature<D>(); break;
			case QUALITY:    current_qcoord_ = asUnsignedInt_(attributes.value("dim")); break;
			case POSITION:   current_pcoord_ = asUnsignedInt_(attributes.value("dim")); break;
			case CONVEXHULL: current_chull_  = new ConvexHullType(); break;
			case HULLPOINT:  hull_position_  = new DPosition<D>(); break;
			case HPOSITION:  current_hcoord_ = asUnsignedInt_(attributes.value("dim")); break;
			case FEATMODEL:
				model_desc_ = new ModelDescription<D>();
		  	param_ = new Param();
		  	if (!attributes.value("name").isEmpty())
		  		model_desc_->setName(attributes.value("name").ascii());
		  	break;
		  case PARAM:
		  	if (!attributes.value("name").isEmpty() && !attributes.value("value").isEmpty() )
		  		param_->setValue(attributes.value("name").ascii(),attributes.value("value").ascii());
		  	break;
		}
		return no_error_;
	}

	template <Size D, typename TraitsT, typename FeatureT>
  bool DFeatureMapHandler<D,TraitsT,FeatureT>::characters( const QString & chars )
  {
		if (is_parser_in_tag_[DESCRIPTION])	// collect Experimental Settings
		{
			exp_sett_ << chars;
			return true;
		}

		// find the tag that the parser is in right now
 		for (Size i=0; i<is_parser_in_tag_.size(); i++)
			if (is_parser_in_tag_[i]){
				switch(i) {
					case FEATINTENSITY: feature_->setIntensity(asDouble_(chars)); break;
					case POSITION:
						feature_->getPosition()[current_pcoord_] = asDouble_(chars);
						break;
					case QUALITY:       feature_->getQuality(current_qcoord_) = asDouble_(chars); break;
					case OVERALLQUALITY:  feature_->getOverallQuality() = asDouble_(chars); break;
					case CHARGE:          feature_->setCharge(asSignedInt_(chars)); break;
					case HPOSITION:       (*hull_position_)[current_hcoord_] = asDouble_(chars); break;
					case META:						feature_->setMetaValue(3,chars.ascii()); break;
				}
			}
		return true;
  }


	template <Size D, typename TraitsT, typename FeatureT>
 	void DFeatureMapHandler<D,TraitsT,FeatureT>::writeTo(std::ostream& os)
	{
		UniqueIdGenerator id_generator = UniqueIdGenerator::instance();

		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
		   << "<featureMap>\n";

		// delegate control to ExperimentalSettings handler
		Internal::MzDataExpSettHandler handler(*((const ExperimentalSettings*)cmap_));
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
			Param::const_iterator piter = modelp.begin();
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
				UnsignedInt hull_size       = current_hull.size();

				for (UnsignedInt j=0;j<hull_size;j++)
				{
					os << "\t\t\t\t<hullpoint>" << std::endl;

					DPosition<D> pos = current_hull[j];
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
