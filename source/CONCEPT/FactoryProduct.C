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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/FactoryProduct.h>

using namespace std;

namespace OpenMS
{
	FactoryProduct::FactoryProduct(): param_(), defaults_(), check_defaults_(true), name_()
	{
	}

	FactoryProduct::FactoryProduct(const FactoryProduct& source)
		: name_(source.name_)
	{
		defaults_ = source.defaults_;
		setParam(source.getParam());
	}

	FactoryProduct::~FactoryProduct()
	{	
	}

	FactoryProduct& FactoryProduct::operator = (const FactoryProduct& source)
	{
		defaults_ = source.defaults_;
		setParam(source.getParam());
		name_ = source.name_;
		return *this;
	}

	FactoryProduct::FactoryProduct(const Param& p)
	{
		FactoryProduct();
		setParam(p);
	}

	void FactoryProduct::setParam(const Param& p)
	{
		if (check_defaults_)
		{
			//cout << "FactoryProduct '" << name_ << "' number of defaults: " << defaults_.size() << endl;
			if (defaults_.size()==0)
			{
				cout << "Warning no default parameters for FactoryProduct '" << name_ << "' specified!" << endl;
			}
			param_ = p;
			param_.setDefaults(defaults_,"",false);
			param_.checkDefaults(defaults_,"");
		}
		else
		{
			param_ = p;
			param_.setDefaults(defaults_);
		}
	}

  const Param& FactoryProduct::getParam() const
	{
		return param_;
	}

  Param& FactoryProduct::getParam()
	{
		return param_;
	}

	const String& FactoryProduct::getName() const
	{
		return name_;
	}
	
	bool FactoryProduct::operator == (const FactoryProduct& rhs) const
	{
			return getParam()==rhs.getParam() && getName()==rhs.getName();
	}

	bool FactoryProduct::operator != (const FactoryProduct& rhs) const
	{
		return !(operator == (rhs));
	}

	std::ostream& operator << (std::ostream& os, const FactoryProduct& prod)
	{
		os << prod.getName() << ":" << std::endl << prod.getParam();
		return os;
	}

}
