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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_XMLSCHEMES_H
#define OPENMS_FORMAT_HANDLERS_XMLSCHEMES_H

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
	/**
		@brief the XML schema namespace

		This namespace contains the XML schemas.
		Schemata are sorted with respect to their version in descending order

		A schema consists of an array of Strings.

		The number of elements has to
		be consistent with the number of Maps in the Handler class.

		The first element has to contain the String the schema is
		recognized by e.g. the name of the xsd-file or the version number.

		All other Strings in the array are used to fill maps in the corresponding
		Handler and contain strings separated by semicolons.
		These strings correspond to values of enumerations (e.g. defined in the PSIOM classes).
		Empty strings are uses if the schema does not define one (e.g. for NULL).
		
		@deprecated This namespace will be removed as soon as SchemaHandler is replaced by XMLHandler (Marc)
	*/
	namespace Schemes
	{
		/// Number of available MzXML schemata
		extern const UInt MzXML_num;

		///Schemata for MzXML
		extern const String MzXML[][11];

		/// Number of available MzData schemata
		extern const UInt MzData_num;

		///Schemata for MzData
		extern const String MzData[][10];

		/// Number of available ExperimentalSettings schemata for MzData
		extern const UInt MzDataExpSett_num;

		///Schemata for the ExperimentalSettings of MzData
		extern const String MzDataExpSett[][18];

		/// Number of available FeatureMap schemata
		extern const UInt FeatureMap_num;

		///Schemata for the FeatureMap
		extern const String FeatureMap[][3];

		/// Number of available DFeaturePairs schemata
		extern const UInt DFeaturePairs_num;

		///Schemata for the FeatureMap
		extern const String DFeaturePairs[][3];

	}
}

#endif // OPENMS_FORMAT_HANDLERS_XMLSCHEMES_H

