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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
	DefaultParamHandler::DefaultParamHandler(const String& name)
		: param_(),
			defaults_(),
			subsections_(),
			error_name_(name),
			check_defaults_(true),
      warn_empty_defaults_(true)
	{
		
	}

	DefaultParamHandler::DefaultParamHandler(const DefaultParamHandler& rhs)
		: param_(rhs.param_),
			defaults_(rhs.defaults_),
			subsections_(rhs.subsections_),
			error_name_(rhs.error_name_),
			check_defaults_(rhs.check_defaults_),
      warn_empty_defaults_(rhs.warn_empty_defaults_)
	{
	}
	
	DefaultParamHandler& DefaultParamHandler::operator= (const DefaultParamHandler& rhs)
	{
		if (&rhs==this) return *this;
		
		//copy members
		param_ = rhs.param_;
		defaults_ = rhs.defaults_;
		subsections_ = rhs.subsections_;
		error_name_ = rhs.error_name_;
		check_defaults_ = rhs.check_defaults_;
    warn_empty_defaults_ = rhs.warn_empty_defaults_;
		
		return *this;
	}
	
	bool DefaultParamHandler::operator== (const DefaultParamHandler& rhs) const
	{
		return 
			param_ == rhs.param_ &&
			defaults_ == rhs.defaults_ &&
			subsections_ == rhs.subsections_ &&
			error_name_ == rhs.error_name_ &&
			check_defaults_ == rhs.check_defaults_ &&
      warn_empty_defaults_ == rhs.warn_empty_defaults_
			;
	}

	DefaultParamHandler::~DefaultParamHandler()
	{
	}	

	void DefaultParamHandler::setParameters(const Param& param)
	{
		//set defaults and apply new parameters
		Param tmp(param);
		tmp.setDefaults(defaults_);
		param_ = tmp;
		
		if (check_defaults_)
		{
			if (defaults_.empty() && warn_empty_defaults_)
			{
				LOG_WARN << "Warning: No default parameters for DefaultParameterHandler '" << error_name_ << "' specified!" << endl;
			}
			
			//remove registered subsections
			for(vector<String>::const_iterator it = subsections_.begin(); it != subsections_.end(); ++it)
			{
				tmp.removeAll(*it+':');
			}
			
			//check defaults
			tmp.checkDefaults(error_name_,defaults_);
		}
		
		//do necessary changes to other member variables
		updateMembers_();	
	}

	void DefaultParamHandler::defaultsToParam_()
	{
		//check if a description is given for all defaults
		bool description_missing = false;
		String missing_parameters;
		for(Param::ParamIterator it = defaults_.begin(); it != defaults_.end();++it)
		{
			//cout << "Name: " << it->getName() << endl;
			if (it->description=="")
			{
				description_missing = true;
				missing_parameters += it.getName()+",";
				break;
			}
		}
		if (description_missing)
		{
			cerr << "Warning: no default parameter description for parameters '" << missing_parameters<< "' of DefaultParameterHandler '" << error_name_ << "' given!" << endl;
		}
		param_.setDefaults(defaults_);
		updateMembers_();	
	}

	void DefaultParamHandler::updateMembers_()
	{
		
	}

	const Param& DefaultParamHandler::getParameters() const
	{
		return param_;
	}
	
	const Param& DefaultParamHandler::getDefaults() const
	{
		return defaults_;
	}

	const String& DefaultParamHandler::getName() const
	{
		return error_name_;
	}

	void DefaultParamHandler::setName(const String& name)
	{
		error_name_ = name;
	}

	const std::vector<String>& DefaultParamHandler::getSubsections() const
	{
		return subsections_;
	}

} // namespace OpenMS


