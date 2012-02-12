// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/CONCEPT/Factory.h>

#include <sstream>

namespace OpenMS
{
	/** 
		@brief Stores the name and parameters of a model.
		
		This class also allows reconstruction of the model.
	
		@see BaseModel
	*/  
  template <UInt D>
	class ModelDescription
	{
		public:
		
		/// Default constructor 
		ModelDescription()
			: name_(), 
				parameters_()
		{
		}
		
		/// copy constructor 
		ModelDescription(const ModelDescription& source)
			:	name_(source.name_),
				parameters_(source.parameters_)
		{
		}

		/// constructor provided for convenience
		ModelDescription(const BaseModel<D>* model)
			: name_(model->getName()),
				parameters_(model->getParameters())
		{
		}

		/// destructor 
		virtual ~ModelDescription()
		{
		}
		
		/// assignment operator 
		virtual ModelDescription& operator = (const ModelDescription& source)
		{
			if (&source == this) return *this;
			
			name_ = source.name_;
			parameters_ = source.parameters_;
			
			return *this;
		}
		
		/// creates model from the parameters defined in this class
		/// returns 0 if no description is set. 
		BaseModel<D>* createModel()
		{
			if (name_ == "") return 0;				
			BaseModel<D>* model = Factory< BaseModel<D> >::create(name_);
			model->setParameters(parameters_);
			return model;
		}
		
		/**	Accessors	*/
		//@{
		/// Non-mutable access to model name
		const String& getName() const
		{ 
			return name_; 
		}
		/// Mutable access to the model name
		String& getName() 
		{ 
			return name_; 
		}
		/// Set the model name
		void setName(const String& name)
		{ 
			name_ = name; 
		}

		/// Non-mutable access to model parameters
		const Param& getParam() const 
		{ 
			return parameters_; 
		}
		/// Mutable access to the model parameters
		Param& getParam() 
		{ 
			return parameters_;
		}
		/// Set the model parameters
		void setParam(const Param& param)
		{ 
			parameters_ = param; 
		}
    
		/**	@name Predicates */
		//@{
		virtual bool operator == (const ModelDescription& rhs) const
		{
			return (name_ == rhs.name_) && (parameters_ == rhs.parameters_);
		}
		
		virtual bool operator != (const ModelDescription& rhs) const
		{
			return !(operator == (rhs));
		}	
		//@}
		
    protected:
    	
    String name_;
		Param parameters_;
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MODELDESCRIPTION_H
