// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>
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
		PTMXMLHandler::PTMXMLHandler( map< String, pair< String, String > >& ptm_informations, const String& filename)
			: XMLHandler(filename,""),
				ptm_informations_(ptm_informations)
		{
		}

		PTMXMLHandler::~PTMXMLHandler()
		{
		}

		void PTMXMLHandler::writeTo(std::ostream& os)
		{
			os << "<PTMs>" << "\n";
			for ( map< String, pair< String, String > >::const_iterator ptm_i = ptm_informations_.begin(); ptm_i != ptm_informations_.end(); ++ptm_i )
			{
					os << "\t<PTM>" << "\n";
					os << "\t\t<name>" << ptm_i->first << "</name>" << "\n"; // see header
					os << "\t\t<composition>" << ptm_i->second.first << "</composition>" << "\n";
					os << "\t\t<possible_amino_acids>" << ptm_i->second.second << "</possible_amino_acids>" << "\n";
					os << "\t</PTM>" << "\n";
			}
			os << "</PTMs>" << "\n";
		}

		void PTMXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& /*attributes*/)
		{
			tag_ = String(sm_.convert(qname)).trim();
			open_tag_ = true;
		}

		void PTMXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const /*qname*/)
		{
// 			tag_ = String(sm_.convert(qname)).trim();
			tag_ = "";
			open_tag_ = false;
		}

		void PTMXMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
		{
			if ( open_tag_ )
			{
				if (tag_ == "name")
				{
					name_ = String(sm_.convert(chars)).trim();
				}
				else if (tag_ == "composition")
				{
					composition_ = String(sm_.convert(chars)).trim();
				}
				else if (tag_ == "possible_amino_acids")
				{
					ptm_informations_[name_] = make_pair(composition_, String(sm_.convert(chars)).trim());
				}
			}
		}

	} // namespace Internal
} // namespace OpenMS
