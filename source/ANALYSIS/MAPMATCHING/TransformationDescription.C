// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

namespace OpenMS
{
	
	TransformationDescription::TransformationDescription()
		: name_(),
		  param_(),
		  pairs_(),
		  trafo_(0)
	{
	}
		
  TransformationDescription::~TransformationDescription()
  {
  	delete trafo_;
  }

	TransformationDescription::TransformationDescription(const TransformationDescription& rhs)
		: name_(rhs.name_),
			param_(rhs.param_),
		  pairs_(rhs.pairs_),
			trafo_(0)
	{
	}

  TransformationDescription& TransformationDescription::operator = (const TransformationDescription& rhs)
  {
    if (this==&rhs) return *this;
    
		name_ = rhs.name_;
		param_ = rhs.param_;
		pairs_ = rhs.pairs_;
		trafo_ = 0;
		
    return *this;
  }

	void TransformationDescription::clear()
	{
		name_ = "";
		param_.clear();
		pairs_.clear();
		delete trafo_;
		trafo_ = 0;
	}

	std::ostream& operator<<(std::ostream& os, TransformationDescription const & td)
	{
		return os <<
		" -- TransformationDescription  BEGIN --\n"
		"name: " << td.getName() << "\n"
		"parameters: " << td.getParameters() <<
		" -- TransformationDescription END --" <<
		std::endl;
	}

} // end of namespace OpenMS

