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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <numeric>

namespace OpenMS
{
    GaussModel::GaussModel()
			: InterpolationModel(), 
				statistics_()
		{
			setName(getProductName());
			
			defaults_.setValue("bounding_box:min",0.0);
			defaults_.setValue("bounding_box:max",1.0);
			defaults_.setValue("statistics:mean",0.0);
			defaults_.setValue("statistics:variance",1.0);
			
			defaultsToParam_();
		}

  	GaussModel::GaussModel(const GaussModel& source)
		: InterpolationModel(source)
		{
			setParam(source.statistics_, source.min_, source.max_);
			updateMembers_();
		}

    GaussModel::~GaussModel()
    {
    }

   	GaussModel& GaussModel::operator = (const GaussModel& source)
		{
			if (&source==this) return *this;
			
			setParam(source.statistics_, source.min_, source.max_);
			InterpolationModel::operator = (source);
			updateMembers_();
			
			return *this;
		}

		void GaussModel::setSamples()
		{
			ContainerType& data = interpolation_.getData();
			data.clear();
			if (max_==min_) return;
			data.reserve( Size ( (max_-min_) / interpolation_step_ + 1 ) );
			CoordinateType pos = min_;
			for ( Size i = 0; pos< max_; ++i)	// old Version: pos<=max_ -> dangerous '==' with double
			{
				pos = min_ + i * interpolation_step_;
				data.push_back( statistics_.normalDensity_sqrt2pi(pos) );
			}
			// scale data so that integral over distribution equals one
			// multiply sum by interpolation_step_ -> rectangular approximation of integral
			IntensityType factor = scaling_ / interpolation_step_ /
										std::accumulate ( data.begin(), data.end(), IntensityType(0) );

			for (ContainerType::iterator it=data.begin();	it!=data.end();	++it) *it *= factor;
			interpolation_.setScale  ( interpolation_step_ );
			interpolation_.setOffset ( min_ );
		}

		void GaussModel::setParam(const BasicStatistics& statistics, CoordinateType min, CoordinateType max)
		{
			statistics_.setMean(statistics.mean());
			statistics_.setVariance(statistics.variance());
			min_ = min;
			max_ = max;

			param_.setValue("bounding_box:min", min_);
			param_.setValue("bounding_box:max", max_);
			param_.setValue("statistics:mean", statistics_.mean());
			param_.setValue("statistics:variance", statistics_.variance());

			setSamples();
		}
		
		void GaussModel::updateMembers_()
		{
			InterpolationModel::updateMembers_();
			
			min_ = param_.getValue("bounding_box:min");
			max_ = param_.getValue("bounding_box:max");
			statistics_.setMean( param_.getValue("statistics:mean") );
			statistics_.setVariance(param_.getValue("statistics:variance"));

			setSamples();
		}

		void GaussModel::setOffset(CoordinateType offset)
		{
			double diff = offset - getInterpolation().getOffset();
			min_ += diff;
			max_ += diff;
			statistics_.setMean(statistics_.mean()+diff);
			
			InterpolationModel::setOffset(offset);

			param_.setValue("bounding_box:min", min_);
			param_.setValue("bounding_box:max", max_);
			param_.setValue("statistics:mean", statistics_.mean());
		}

		const GaussModel::CoordinateType GaussModel::getCenter() const
		{
			return statistics_.mean();
		}

}
