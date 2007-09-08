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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_FEATUREPAIRSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_FEATUREPAIRSHANDLER_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/FORMAT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>

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
  class FeaturePairsHandler
		: public SchemaHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      ///
      FeaturePairsHandler(std::vector< ElementPair < Feature > > & map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				pairs_(&map), cpairs_(0),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()
		  {
				fillMaps_(Schemes::DFeaturePairs[schema_]);
				setMaps_(TAGMAP, ATTMAP);
			}

      ///
      FeaturePairsHandler(const std::vector< ElementPair < Feature > > & map, const String& filename)
      : SchemaHandler(TAG_NUM,MAP_NUM,filename),
				pairs_(0), cpairs_(&map),
				id_generator_(UniqueIdGenerator::instance()),
				pair_(), feature_()
		  { 
				fillMaps_(Schemes::DFeaturePairs[schema_]);
				setMaps_(TAGMAP, ATTMAP);
			}
      ///
      virtual ~FeaturePairsHandler()
      {
      }
      //@}

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
      virtual void characters(const XMLCh* const chars, unsigned int length);

		/// Print the contents to a stream
			void writeTo(std::ostream& os);

		protected:

		/** @brief indices for tags used by FeatureXMLFile

			Used to access is_parser_in_tag_.
			If you add tags, also add them to XMLSchemes.h.
			Add no elements to the enum after TAG_NUM.
		*/
		enum Tags {	TAGNULL, PAIRLIST, PAIR, PAIRQUALITY, FIRST, SECOND,
								FEATURE, POSITION, FEATINTENSITY, QUALITY,
								OVERALLQUALITY, CHARGE, FEATMODEL, PARAM, CONVEXHULL,
								HULLPOINT, HPOSITION, TAG_NUM};

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

		/// Vector of pairs to be read
		std::vector< ElementPair < Feature > > * pairs_;
		/// Vector of pairs to be written
		const std::vector< ElementPair < Feature > > * cpairs_;
		/// ID generator
		UniqueIdGenerator id_generator_;

		/// The current coordinates
		UInt current_pcoord_;
		UInt current_qcoord_;
		UInt current_hcoord_;

		// temporary datastructures to hold parsed data
		ElementPair < Feature >* pair_;
		Feature* feature_;
		ModelDescription<2>* model_desc_;
		Param* param_;
		ConvexHull2D* current_chull_;
		Feature::PositionType* hull_position_;

		void writeFeature_(std::ostream& os, Feature dfeat);

	}; // end of class FeaturePairsHandler

	//--------------------------------------------------------------------------------

 	} // namespace Internal
} // namespace OpenMS

#endif
