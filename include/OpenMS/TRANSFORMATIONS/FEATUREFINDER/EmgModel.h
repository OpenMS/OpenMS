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


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


namespace OpenMS
{
	/** 
		@brief Exponentially modified gaussian distribution model for elution profiles.
    
    @htmlinclude OpenMS_EmgModel.parameters
	*/
	class OPENMS_DLLAPI EmgModel
		: public InterpolationModel
	{

	 public:
		typedef InterpolationModel::CoordinateType CoordinateType;
		typedef Math::BasicStatistics<CoordinateType > BasicStatistics;
		typedef LinearInterpolation::container_type ContainerType;

		/// Default constructor
		EmgModel();

		/// copy constructor
		EmgModel(const EmgModel& source);

		/// destructor
		virtual ~EmgModel();

		/// assignment operator
		virtual EmgModel& operator = (const EmgModel& source);

		/// create new EmgModel object (needed by Factory)
		static BaseModel<1>* create()
		{
			return new EmgModel();
  	}

		/// name of the model (needed by Factory)
		static const String getProductName()
		{
			return "EmgModel";
		}

		/// set offset without being computing all over and without any discrepancy
    void setOffset(CoordinateType offset);

		/// set sample/supporting points of interpolation
		void setSamples();

		/// get the center of the Gaussian model i.e. the position of the maximum
		CoordinateType getCenter() const;

	 protected:
		CoordinateType  min_;
		CoordinateType  max_;
		BasicStatistics statistics_;
		CoordinateType height_;
		CoordinateType width_;
		CoordinateType symmetry_;
		CoordinateType retention_;
		
		void updateMembers_();
	};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EMGMODEL_H
