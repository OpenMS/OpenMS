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

#ifndef OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H
#define OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H

#include <OpenMS/FORMAT/Param.h>

#include <vector>

namespace OpenMS
{
	/**	
		@brief A base class for all classes handling default parameters.
		
		This class facilitates the handling of parameters:
		- it manages default parameter (defaults_)
		- it checks for wrong/misspelled parameters
		- subsections that are i.e. passed to other classes can be excluded from the check (subsections_)
		- it keeps member variables in synronicity with the parameters stored in param_
		
		Extra member variables are needed if getting the value from param_ would be too slow
		e.g. when they are used in methods that are called very often.
		
		No matter if you have extra variables or not, do the following:
		- Set defaults_ and subsections_ in the derived classes' default constructor.
		- Call defaultsToParam_() at the end of derived classes' default constructor.
			It copies the defaults to param_ (and calls updateMembers_()).
		
		If you have extra member variables you need to syncronize with param_, do the following:
		- Implement the updateMembers_() method. It is used after each change of param_
			in order to update the extra member variables. If the base class is a DefaultParamHandler as well
			make sure to call the updateMembers_() method of the base class in the updateMembers_() method.
		- Call updateMembers_() at the end of the derived classes' copy constructor and assignment operator.
		- If you need mutable access to the extra member variables, provide a set-method and make sure to set
		  the corresponding value in param_ as well!

		@ingroup Datastructures
	*/
	class DefaultParamHandler
	{
		public:
			/// Constructor with name that is diplayed in error messages
			DefaultParamHandler(const String& name);
			
			/// Copy constructor
			DefaultParamHandler(const DefaultParamHandler& rhs);

			/// Destructor
			virtual ~DefaultParamHandler();
			
			///@brief Assignment operator.
			virtual DefaultParamHandler& operator= (const DefaultParamHandler& rhs);
			
			/// Equality operator
			virtual bool operator== (const DefaultParamHandler& rhs) const;
			
			/**
				@brief Sets the parameters
				
				It also applies the default parameters to @p param and checks for unknown parameters.
			*/
			void setParameters(const Param& param);
			
			/// Non-mutable access to the parameters
			const Param& getParameters() const;
			
			/// Non-mutable access to the default parameters
			const Param& getDefaults() const;			
			
			/// Non-mutable access to the name
    	const String& getName() const;

			/// Mutable access to the name
    	void setName(const String& name);
    	
    	/// Non-mutable access to the registered subsections
    	const std::vector<String>& getSubsections() const;
			
		protected:
			/**
				@brief This method is used to update extra member variables at the end of the setParam() method.
				
				Also call it at the end of the derived classes' copy constructor and assignment operator.
				
				The default implementation is empty.
			*/
			virtual void updateMembers_();	
			
			///Updates the parameters after the defaults have been set in the constructor
			void defaultsToParam_();
			
			///Container for current parameters
			Param param_;
			
			/**
				@brief Container for default paramters. This member should be filled in the constructor of derived classes!
				
				@note call the setParam_() method at the end of the constructor in order to copy the defaults to param_.
			*/
			Param defaults_;
			
			/**
				@brief Container for registered subsections. This member should be filled in the constructor of derived classes!
				
				@note Do not add a ':' character at the end of subsections.
			*/
			std::vector<String> subsections_;
			
			/// Name that is displayed in error messages during the parameter checking
			String error_name_;
			
			/**
				@brief If this member is set to false no checking if parameters in done;
				
				The only reason to set this member to false is that the derived class has no parameters!
			*/
			bool check_defaults_;

		private:
			/// Hidden as a name is required!
			DefaultParamHandler();
	
	}; //class

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DEFAULTPARAMHANDLER_H
