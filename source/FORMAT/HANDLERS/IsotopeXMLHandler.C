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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/IsotopeXMLHandler.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <xercesc/sax2/Attributes.hpp>

#include <iostream>
#include <map>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
		IsotopeXMLHandler::IsotopeXMLHandler(
			map< String, vector< pair< DoubleReal, DoubleReal > > >& isotope_informations,
			const String& filename):
				XMLHandler(filename),
				isotope_informations_(isotope_informations),
				const_isotope_informations_()
		{}

		IsotopeXMLHandler::IsotopeXMLHandler(
			const map< String, vector< pair< DoubleReal, DoubleReal > > >& isotope_informations,
			const String& filename):
				XMLHandler(filename),
				isotope_informations_(),
				const_isotope_informations_(isotope_informations)
		{}

		IsotopeXMLHandler::~IsotopeXMLHandler()
		{}

		void IsotopeXMLHandler::writeTo(ostream& os)
// 		throw (Exception::IndexOverflow)
		{
			os << "<isotopes>" << endl;
			for ( map< String, vector< pair< DoubleReal, DoubleReal > > >::const_iterator element_i = const_isotope_informations_.begin(); element_i != const_isotope_informations_.end(); ++element_i )
			{
				os << "\t<element>" << endl;
				os << "\t\t<symbol>" << element_i->first << "</symbol>" << endl;
				for ( vector< pair< DoubleReal, DoubleReal > >::const_iterator isotope_i = element_i->second.begin(); isotope_i != element_i->second.end(); ++isotope_i )
				{
					os << "\t\t<isotope>" << endl;
					os << "\t\t\t<mass>" << isotope_i->first << "</mass>" << endl;
					os << "\t\t\t<probability>" << isotope_i->second << "</probability>" << endl;
					os << "\t\t</isotope>" << endl;
				}
				os << "\t</element>" << endl;
			}
			os << "</isotopes>" << endl;
		}

		void IsotopeXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& /*attributes*/)
		{
			tag_ = String(xercesc::XMLString::transcode(qname)).trim();
			open_tag_ = true;
		}

		void IsotopeXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const /*qname*/)
		{
// 			tag_ = String(xercesc::XMLString::transcode(qname)).trim();
			tag_ = "";
			open_tag_ = false;
		}

		void IsotopeXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
		{
			if ( open_tag_ )
			{
				if (tag_ == "symbol")
				{
					symbol_ = String(xercesc::XMLString::transcode(chars)).trim();
					isotope_informations_[symbol_] = vector< pair< DoubleReal, DoubleReal > >();
				}
				else if (tag_ == "mass")
				{
					mass_ = String(xercesc::XMLString::transcode(chars)).toDouble();
				}
				else if (tag_ == "probability")
				{
					isotope_informations_[symbol_].push_back(make_pair(mass_, String(xercesc::XMLString::transcode(chars)).toDouble()));
				}
			}
		}

	} // namespace Internal
} // namespace OpenMS
