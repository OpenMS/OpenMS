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
// $Id: MetaInfoRegistry.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

class MetaInfoRegistry
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/MetaInfoRegistry.h>
	using namespace OpenMS;
%End

  public:
    
	class UnregisteredName	
		: Exception::Base
	{
		public:
		UnregisteredName(const char* /file/, int /line/, const std::string& /s/) throw();
		~UnregisteredName() throw();
	};

	class UnregisteredIndex	
		: Exception::Base
	{
		public:
		UnregisteredIndex(const char* /file/, SignedInt /line/, SignedInt /index/) throw();
		~UnregisteredIndex() throw();
	};
    
	MetaInfoRegistry();
	MetaInfoRegistry(const MetaInfoRegistry& /rhs/);
	~MetaInfoRegistry();
	// MetaInfoRegistry& operator = (const MetaInfoRegistry& rhs);
	UnsignedInt registerName(const std::string& /name/, const std::string& /description/, const std::string& /unit/ = "");
	UnsignedInt getIndex(const std::string& /name/) const;
	std::string getDescription(UnsignedInt /index/) const throw(MetaInfoRegistry::UnregisteredIndex);
	std::string getDescription(const std::string& /name/) const throw(MetaInfoRegistry::UnregisteredName);
	std::string getUnit(UnsignedInt /index/) const throw(MetaInfoRegistry::UnregisteredIndex);
	std::string getUnit(const std::string& /name/) const throw(MetaInfoRegistry::UnregisteredName);
};

