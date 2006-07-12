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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_DFEATUREPAIRSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_DFEATUREPAIRSHANDLER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DPosition.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>

// STL includes
#include <iostream>
#include <valarray>
#include <string>

namespace OpenMS
{
	namespace Internal
	{

	/** @brief XML Handler for a DFeaturePairVector i.e. a vector consisting
	    of pairs of features as there are produced by the FeatureMatcher.	 
	 */
  template <Size D, typename FeatureT = DFeature<D> >
  class DFeaturePairsHandler
		: public SchemaHandler
  {
    public:
    	/**
				@name Type definitions
			*/
			//@{
			typedef FeatureT FeatureType;
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
			//@}
						
      /**@name Constructors and destructor */
      //@{
      ///
      DFeaturePairsHandler(DFeaturePairVector<D,FeatureType>& map)
      : SchemaHandler(TAG_NUM,MAP_NUM),
				pairs_(&map), cpairs_(0),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()		
  		{
				file_ = __FILE__;
				fillMaps_(Schemes::DFeaturePairs[schema_]);
			}
      
      ///
      DFeaturePairsHandler(const DFeaturePairVector<D,FeatureType>& map)
      : SchemaHandler(TAG_NUM,MAP_NUM),
				pairs_(0), cpairs_(&map),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()	
  		{
				file_ = __FILE__;
				fillMaps_(Schemes::DFeaturePairs[schema_]);
			}
      ///
      virtual ~DFeaturePairsHandler() 
      {
      }
      //@}

      /// This function is called for each closing tag in the XML file.
      virtual bool endElement( const QString & uri, const QString & local_name,
															 const QString & qname );

      /// This function is called for each opening XML tag in the file.
      virtual bool startElement(const QString & uri, const QString & local_name,
																const QString & qname, const QXmlAttributes & attributes );

			/// This function is called by the parser for each chunk of
		  /// characters between two tags.
      virtual bool characters( const QString & chars );

  		/// Print the contents to a stream
			void writeTo(std::ostream& os);

		protected:

		/** @brief indices for tags used by DFeatureMapFile

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags {	TAGNULL, PAIRLIST, PAIR, PAIRQUALITY, FIRST, SECOND,
								FEATURE, POSITION, FEATINTENSITY, QUALITY,
								OVERALLQUALITY, CHARGE, FEATMODEL, PARAM, CONVEXHULL,
								HULLPOINT, HPOSITION, TAG_NUM};

		/** @brief indices for enum2str-maps used by DFeatureMapFile

			Used to access enum2str_().
			If you add maps, also add them to XMLSchemes.h.
			Add no elements to the enum after MAP_NUM.
			Each map corresponds to a string in XMLSchemes.h.
		*/
		enum MapTypes {	TAGMAP, MAP_NUM };

		/// Vector of pairs to be read
		DFeaturePairVector<D,FeatureType>* pairs_;
		/// Vector of pairs to be written
		const DFeaturePairVector<D,FeatureType>* cpairs_;
		/// ID generator
		UniqueIdGenerator id_generator_; 

		/// The current coordinates
		UnsignedInt current_pcoord_;
		UnsignedInt current_qcoord_;
		UnsignedInt current_hcoord_;
			
		// temporary datastructures to hold parsed data
		DFeaturePair<D, FeatureType>* pair_;
		FeatureType* feature_;
		ModelDescription<D>* model_desc_;
		Param* param_;
		ConvexHullType* current_chull_;
		DPosition<D>* hull_position_;

		void writeFeature_(std::ostream& os, FeatureType dfeat);

	}; // end of class DFeaturePairsHandler

	//--------------------------------------------------------------------------------

	template <Size D, typename FeatureT>
  bool DFeaturePairsHandler<D,FeatureT>::startElement(
				const QString & /*uri*/, const QString & /*local_name*/,
				const QString & qname, const QXmlAttributes & attributes)
	{
		int tag = str2enum_(TAGMAP,qname,"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;

		switch(tag)
		{
			case FEATURE: 	 feature_        = new DFeature<D>(); break;
			case PAIR:		 	 pair_	         = new DFeaturePair<D>(); break;
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
		return true;
	}

  template <Size D, typename FeatureT>
  bool DFeaturePairsHandler<D,FeatureT>::endElement( const QString & /*uri*/,
	const QString & /*local_name*/,	 const QString & qname )
 	{
 		int tag = str2enum_(TAGMAP,qname,"closing tag");  // index of current tag
		is_parser_in_tag_[tag] = false;
		switch(tag) {
			case FIRST:
				pair_->setFirst(*feature_);
				delete feature_;
				break;
			case SECOND:
				pair_->setSecond(*feature_);
				delete feature_;
				break;
			case PAIR:
				pairs_->push_back(*pair_);
				delete pair_;
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

	template <Size D, typename FeatureT>
  bool DFeaturePairsHandler<D,FeatureT>::characters( const QString & chars )
  {
		for (Size i=0; i<is_parser_in_tag_.size(); i++)
		{
			if (is_parser_in_tag_[i])
			{
				switch(i)
				{
					case FEATINTENSITY:		feature_->setIntensity(asDouble_(chars)); break;
					case POSITION:        feature_->getPosition()[current_pcoord_] = asDouble_(chars); break;
					case QUALITY:         feature_->getQuality(current_qcoord_) = asDouble_(chars); break;
					case OVERALLQUALITY:  feature_->getOverallQuality() = asDouble_(chars); break;
					case CHARGE:          feature_->setCharge(asSignedInt_(chars)); break;
					case HPOSITION:       (*hull_position_)[current_hcoord_] = asDouble_(chars); break;
					case PAIRQUALITY:			pair_->setQuality(asDouble_(chars));
				}
			}
   	}
		return true;
  }

	template <Size D, typename FeatureT>
	void DFeaturePairsHandler<D,FeatureT>::writeTo(std::ostream& os)
	{

		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
		os << "<pairlist>" << std::endl;

		// write features with their attributes
		for (UnsignedInt s=0; s<cpairs_->size(); s++)
		{
			const DFeaturePair<D>& pair = (*cpairs_)[s];

			os << "<pair nr=\"" << s << "\">" << std::endl;
			os << "\t<pairquality>" << pair.getQuality() << "</pairquality>" << std::endl;

			os << "\t<first>" << std::endl;
			FeatureType first = pair.getFirst();
			writeFeature_(os,first);
			os << "\t</first>" << std::endl;

			os << "\t<second>" << std::endl;
			FeatureType seco  = pair.getSecond();
			writeFeature_(os,seco);
			os << "\t</second>" << std::endl;

			os << "</pair>" << std::endl;

		} // end for ( features )

		os << "</pairlist>" << std::endl;
		os <<
			"<!-- Local Variables: -->\n"
			"<!-- mode: nxml -->\n"
			"<!-- tab-width: 2 -->\n"
			"<!-- End: -->\n"
			;
	}



	template <Size D, typename FeatureT>
 	void DFeaturePairsHandler<D,FeatureT>::writeFeature_(std::ostream& os, FeatureType dfeat)
 	{
		os << "\t<feature id=\"" << id_generator_.getUID() << "\">" << std::endl;

		DPosition<D> pos = dfeat.getPosition();
		UnsignedInt dpos_size = pos.size();

		for (UnsignedInt i=0; i<dpos_size;i++)
		{
			os <<	"\t\t<position dim=\"" << i << "\">" << pos[i] << "</position>" << 	std::endl;
		}

		os << "\t\t<intensity>" << dfeat.getIntensity() << "</intensity>" << std::endl;

		for (UnsignedInt i=0; i<dpos_size;i++)
		{
			os << "\t\t<quality dim=\"" << i << "\">" << dfeat.getQuality(i) << "</quality>" << std:: endl;
		}

		os << "\t\t<overallquality>" << dfeat.getOverallQuality() << "</overallquality>" << std:: endl;
		os << "\t\t<charge>" << dfeat.getCharge() << "</charge>" << std:: endl;

		// write model description
		ModelDescription<D> desc = dfeat.getModelDescription();
		os << "\t\t<model name=\"" << desc.getName() << "\">" << std:: endl;
		Param modelp = desc.getParam();
		Param::ConstIterator piter = modelp.begin();
		while (piter != modelp.end())
		{
			os << "\t\t\t<param name=\"" << piter->first << "\" value=\"" << piter->second << "\">";
			os << "</param>" << std::endl;
			piter++;
		}
		os << "\t\t</model>" << std::endl;

		// write convex hull
		DFeature<2>::ConvexHullVector hulls = dfeat.getConvexHulls();
		DFeature<2>::ConvexHullVector::iterator citer = hulls.begin();

		UnsignedInt hulls_count = hulls.size();

		for (UnsignedInt i=0;i<hulls_count; i++)
		{
			os << "\t\t<convexhull nr=\"" << i << "\">" << std:: endl;

			ConvexHullType current_hull = hulls[i];
			UnsignedInt hull_size       = current_hull.size();

			for (UnsignedInt j=0;j<hull_size;j++)
			{
				os << "\t\t\t<hullpoint>" << std::endl;

				DPosition<D> pos = current_hull[j];
				UnsignedInt pos_size = pos.size();
				for (UnsignedInt k=0; k<pos_size; k++)
				{
					os << "\t\t\t\t<hposition dim=\"" << k << "\">" << pos[k] << "</hposition>" << std::endl;
				}

				os << "\t\t\t</hullpoint>" << std::endl;
			} // end for (..hull_size..)

			os << "\t\t</convexhull>" << std::endl;
		} // end  for ( ... hull_count..)

		os << "\t</feature>\n";
 	}


	} // namespace Internal
} // namespace OpenMS

#endif
