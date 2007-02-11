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

#ifndef OPENMS_FORMAT_HANDLERS_DFEATUREPAIRSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_DFEATUREPAIRSHANDLER_H

#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
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

#include <xercesc/sax2/Attributes.hpp>

namespace OpenMS
{
	namespace Internal
	{

	/** @brief XML Handler for a DFeaturePairVector
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
      DFeaturePairsHandler(DFeaturePairVector<D,FeatureType>& map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				pairs_(&map), cpairs_(0),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()
		{
				fillMaps_(Schemes::DFeaturePairs[schema_]);
				setMaps_(TAGMAP, ATTMAP);
			}

      ///
      DFeaturePairsHandler(const DFeaturePairVector<D,FeatureType>& map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				pairs_(0), cpairs_(&map),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()
		{
				fillMaps_(Schemes::DFeaturePairs[schema_]);
				setMaps_(TAGMAP, ATTMAP);
			}
      ///
      virtual ~DFeaturePairsHandler()
      {
      }
      //@}

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
      virtual void characters(const XMLCh* const chars, const unsigned int length);

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
  void DFeaturePairsHandler<D,FeatureT>::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		int tag = enterTag(qname, attributes);

		String tmp_str;
		switch(tag)
		{
			case FEATURE:	 feature_        = new DFeature<D>(); break;
			case PAIR:			 pair_	         = new DFeaturePair<D>(); break;
			case QUALITY:
				tmp_str = getAttributeAsString_(DIM);
				current_qcoord_ = asUnsignedInt_(tmp_str);
				break;
			case POSITION:
				tmp_str = getAttributeAsString_(DIM);
				current_pcoord_ = asUnsignedInt_(tmp_str);
				break;
		case CONVEXHULL: current_chull_  = new ConvexHullType(); break;
		case HULLPOINT:  hull_position_  = new DPosition<D>(); break;
		case HPOSITION:
				tmp_str = getAttributeAsString_(DIM);
			current_hcoord_ = asUnsignedInt_(tmp_str);
			break;
		case FEATMODEL:
			model_desc_ = new ModelDescription<D>();
			param_ = new Param();
			tmp_str = getAttributeAsString_(NAME);
			if (tmp_str != "")
				model_desc_->setName(tmp_str);
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

  template <Size D, typename FeatureT>
  void DFeaturePairsHandler<D,FeatureT>::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
	  int tag = leaveTag(qname);
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
				current_chull_->addPoint(*hull_position_);
				delete hull_position_;
				break;
			case CONVEXHULL:
				feature_->getConvexHulls().push_back(*current_chull_);
				delete current_chull_;
				break;
		}
  }

	template <Size D, typename FeatureT>
  void DFeaturePairsHandler<D,FeatureT>::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
		for (Size i=0; i<is_parser_in_tag_.size(); i++)
		{
			if (is_parser_in_tag_[i])
			{
				switch(i)
				{
					case FEATINTENSITY:		feature_->setIntensity(asDouble_(xercesc::XMLString::transcode(chars))); break;
					case POSITION:        feature_->getPosition()[current_pcoord_] = asDouble_(xercesc::XMLString::transcode(chars)); break;
					case QUALITY:         feature_->getQuality(current_qcoord_) = asDouble_(xercesc::XMLString::transcode(chars)); break;
					case OVERALLQUALITY:  feature_->getOverallQuality() = asDouble_(xercesc::XMLString::transcode(chars)); break;
					case CHARGE:          feature_->setCharge(asSignedInt_(xercesc::XMLString::transcode(chars))); break;
					case HPOSITION:       (*hull_position_)[current_hcoord_] = asDouble_(xercesc::XMLString::transcode(chars)); break;
					case PAIRQUALITY:			pair_->setQuality(asDouble_(xercesc::XMLString::transcode(chars)));
				}
			}
	}
  }

	template <Size D, typename FeatureT>
	void DFeaturePairsHandler<D,FeatureT>::writeTo(std::ostream& os)
	{

		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
		os << "<featurePairs>" << std::endl;

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

		os << "</featurePairs>" << std::endl;
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
			os <<	"\t\t<position dim=\"" << i << "\">" << pos[i] << "</position>" <<	std::endl;
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
			UnsignedInt hull_size = current_hull.getPoints().size();

			for (UnsignedInt j=0;j<hull_size;j++)
			{
				os << "\t\t\t<hullpoint>" << std::endl;

				DPosition<D> pos = current_hull.getPoints()[j];
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


// Please leave the page-feed character (Ctrl-L) in the following line,
// otherwise emacs will not select c++ mode when opening this file.

// Thanks!

// EOF
