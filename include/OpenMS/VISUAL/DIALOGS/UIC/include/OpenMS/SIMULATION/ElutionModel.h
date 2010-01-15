// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_ELUTIONMODEL_H
#define OPENMS_SIMULATION_ELUTIONMODEL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

#include <boost/math/tr1.hpp>

namespace OpenMS
{
	/**
		@brief Exponentially modified gaussian distribution model for elution profiles.

    @htmlinclude OpenMS_ElutionModel.parameters

	*/
	class OPENMS_DLLAPI ElutionModel
		: public InterpolationModel
	{

	 public:
		typedef InterpolationModel::CoordinateType CoordinateType;
		typedef Math::BasicStatistics<CoordinateType > BasicStatistics;
    typedef LinearInterpolation::container_type ContainerType;

		/// Default constructor
		ElutionModel();

		/// copy constructor
		ElutionModel(const ElutionModel& source);

		/// destructor
		virtual ~ElutionModel();

		/// assignment operator
		virtual ElutionModel& operator = (const ElutionModel& source);

		/// create new ElutionModel object (needed by Factory)
		static BaseModel<1>* create()
		{
			return new ElutionModel();
  	}

		/// name of the model (needed by Factory)
		static const String getProductName()
		{
			return "ElutionModel";
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

} // namespace OpenMS

#endif // OPENMS_SIMULATION_ELUTIONMODEL_H
