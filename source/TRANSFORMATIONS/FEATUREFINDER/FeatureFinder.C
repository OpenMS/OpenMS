// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModelFitter.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	FeatureFinder::FeatureFinder()
	{
		traits_ = new FeaFiTraits();
	}
	
	FeatureFinder::FeatureFinder(const FeatureFinder& source)
		: traits_(source.traits_),
		seeders_(source.seeders_),
		extenders_(source.extenders_),
		fitters_(source.fitters_)
	{
	}
	
	FeatureFinder::~FeatureFinder()
	{
	}
	
	FeatureFinder& FeatureFinder::operator = (const FeatureFinder& source)
	{
		if (&source == this) return *this;
	
		traits_ = source.traits_;
		seeders_ = source.seeders_;
		extenders_ = source.extenders_;
		fitters_ = source.fitters_;
	
		return *this;
	}
	
	bool FeatureFinder::operator == (const FeatureFinder& rhs) const
	{
		return (traits_	   == rhs.traits_ &&
				   seeders_   == rhs.seeders_ &&
				   extenders_ == rhs.extenders_ &&
				   fitters_	   == rhs.fitters_);
	}
	
	void FeatureFinder::addSeeder(const String& name, const Param* param)
	{
		BaseSeeder* p = Factory<BaseSeeder>::create(name);
		if (param) p->setParameters(*param);
		seeders_.push_back(p);
	}
	
	void FeatureFinder::addExtender(const String& name, const Param* param)
	{
		BaseExtender* p = Factory<BaseExtender>::create(name);
		if (param) p->setParameters(*param);
		extenders_.push_back(p);
	}
	
	void FeatureFinder::addFitter(const String& name, const Param* param)
	{
		BaseModelFitter* p = Factory<BaseModelFitter>::create(name);
		if (param) p->setParameters(*param);
		fitters_.push_back(p);
	}
	
	void FeatureFinder::removeSeeder(const String& name)
	{
		for (SeederVector::iterator it=seeders_.begin(); it!=seeders_.end(); it++)
		{
			if ((*it)->getName() == name)
			{
				seeders_.erase(it);
				break;
			}
		}
	}
	
	void FeatureFinder::removeExtender(const String& name)
	{
		for (ExtenderVector::iterator it=extenders_.begin(); it!=extenders_.end(); it++)
		{
			if ((*it)->getName() == name)
			{
				extenders_.erase(it);
				break;
			}
		}
	}
	
	void FeatureFinder::removeFitter(const String& name)
	{
		for (FitterVector::iterator it=fitters_.begin(); it!=fitters_.end(); it++)
		{
			if ((*it)->getName() == name)
			{
				fitters_.erase(it);
				break;
			}
		}
	}
	
	bool FeatureFinder::setParam(const Param& param)
	{
		param_ = param;
		if (param_.empty())
		{
			return false;
		}
		else
		{
			if (!setModule("Seeders"))
			{
				cerr << "Error: Could not set Seeders module" << endl;
				return false;
			}
			if (!setModule("Extenders"))
			{
				cerr << "Error: Could not set Extenders module" << endl;
				return false;
			}
			if (!setModule("ModelFitters"))
			{
				cerr << "Error: Could not set ModelFitters module" << endl;
				return false;
			}
		}
		return true;
	}
	
	
	const DFeatureMap<2>& FeatureFinder::run()
	{
		if (!traits_)
		{
			throw Exception::Base(__FILE__, __LINE__,"FeatureFinder::run()","NotInitialized","FeatureFinder has not been initialized");
		}
		
		for (SeederVector::iterator it=seeders_.begin(); it!=seeders_.end(); ++it)
		{
			(*it)->setTraits(traits_);
		}
		for (ExtenderVector::iterator it=extenders_.begin(); it!=extenders_.end(); ++it)
		{
			(*it)->setTraits(traits_);
		}
		for (FitterVector::iterator it=fitters_.begin(); it!=fitters_.end(); ++it)
		{
			(*it)->setTraits(traits_);
		}
		
		return traits_->run(seeders_,extenders_,fitters_);
	}
	
	std::ostream& operator << (std::ostream& os, const FeatureFinder& finder)
	{
		os << "Seeders: ";
	
		for (FeatureFinder::SeederVector::const_iterator it = finder.seeders_.begin(); it != finder.seeders_.end(); ++it)
		{
			os << (*it)->getName() << " ";
		}
		
		os << "\nExtenders: ";
		for (FeatureFinder::ExtenderVector::const_iterator it = finder.extenders_.begin(); it != finder.extenders_.end(); ++it)
		{
			os << (*it)->getName() << " ";
		}
		
		os << "\nModelfitters: ";
		for (FeatureFinder::FitterVector::const_iterator it = finder.fitters_.begin(); it != finder.fitters_.end(); ++it)
		{
			os << (*it)->getName() << " ";
		}
		
		return os;
	}

}
