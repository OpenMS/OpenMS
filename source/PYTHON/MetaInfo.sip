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
// $Id: MetaInfo.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

class MetaInfo
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/MetaInfo.h>
	using namespace OpenMS;
%End
    public:
		MetaInfo();
		MetaInfo(const MetaInfo& /rhs/);
		~MetaInfo();
		// MetaInfo& operator = (const MetaInfo& /rhs/);
		
		/// returns the value corresponding to a  string
		DataValue getValue(const std::string& /name/) const throw (MetaInfoRegistry::UnregisteredName);
		/// returns the value corresponding to an index
		DataValue getValue(UnsignedInt /index/) const throw (MetaInfoRegistry::UnregisteredIndex);
		
		bool exists(const std::string& /name/);
		bool exists(UnsignedInt /index/);
		void setValue(const std::string& /name/, const std::string& /value/);
		void setValue(UnsignedInt /index/, const std::string& /value/);
		void setValue(const std::string& /name/, SignedInt /value/);
		void setValue(UnsignedInt /index/, const SignedInt /value/);
		void setValue(const std::string& /name/, double /value/);
		void setValue(UnsignedInt /index/, double /value/);
		void setValue(const std::string& /name/, const DataValue& /value/);
		void setValue(UnsignedInt /index/, const DataValue& /value/);	
		static MetaInfoRegistry& registry();
};
