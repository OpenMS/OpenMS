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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <numeric>

namespace OpenMS
{
    BiGaussModel::BiGaussModel()
		: InterpolationModel<>(), statistics1_(), statistics2_()
		{
			setName(getProductName());
			
			defaults_.setValue("bounding_box:min",0.0);
			defaults_.setValue("bounding_box:max",1.0);
			defaults_.setValue("statistics:mean",0.0);
			defaults_.setValue("statistics:variance1",1.0);
			defaults_.setValue("statistics:variance2",1.0);
			
			defaultsToParam_();
		}

  	BiGaussModel::BiGaussModel(const BiGaussModel& source)
		: InterpolationModel<>(source)
		{
			setParam(source.statistics1_.mean(), source.statistics1_.variance(),
							 source.statistics2_.variance(),source.min_, source.max_);
			updateMembers_();
		}

    BiGaussModel::~BiGaussModel()
    {
    }

   	BiGaussModel& BiGaussModel::operator = (const BiGaussModel& source)
		{
			if (&source == this) return *this;
			
			InterpolationModel<>::operator = (source);
			setParam(source.statistics1_.mean(), source.statistics1_.variance(),
							 source.statistics2_.variance(),source.min_, source.max_);
			updateMembers_();
			
			return *this;
		}

		void BiGaussModel::setSamples()
		{
			ContainerType& data = interpolation_.getData();
			data.clear();
			if (max_==min_) return;
			data.reserve( Size ( (max_-min_) / interpolation_step_ + 1 ) );
			CoordinateType pos = min_;

			for ( Size i = 0; pos< max_; ++i)
			{
				pos = min_ + i * interpolation_step_;
				if (pos < statistics1_.mean())
					data.push_back( statistics1_.normalDensity_sqrt2pi(pos) );
				else
					data.push_back( statistics2_.normalDensity_sqrt2pi(pos) );
			}
			// scale data so that integral over distribution equals one
			// multiply sum by interpolation_step_ -> rectangular approximation of integral
			IntensityType factor = scaling_ / interpolation_step_ /
						std::accumulate ( data.begin(), data.end(), IntensityType(0) );

			for (ContainerType::iterator it=data.begin();	it!=data.end();	++it){
				*it *= factor;
				}

			interpolation_.setScale  ( interpolation_step_ );
			interpolation_.setOffset ( min_ );
		}

		void BiGaussModel::setParam(CoordinateType mean, CoordinateType variance1, CoordinateType variance2,	CoordinateType min, CoordinateType max)
		{
			min_ = min;
			max_ = max;
			statistics1_.setMean(mean);
			statistics2_.setMean(mean);
			statistics1_.setVariance(variance1);
			statistics2_.setVariance(variance2);

			param_.setValue("bounding_box:min", min_);
			param_.setValue("bounding_box:max", max_);
			param_.setValue("statistics:mean", statistics1_.mean());
			param_.setValue("statistics:variance1", statistics1_.variance());
			param_.setValue("statistics:variance2", statistics2_.variance());

			setSamples();
		}

		void BiGaussModel::updateMembers_()
		{
			InterpolationModel<>::updateMembers_();
			
			min_ = param_.getValue("bounding_box:min");
			max_ = param_.getValue("bounding_box:max");
			statistics1_.setMean(param_.getValue("statistics:mean"));
			statistics2_.setMean(param_.getValue("statistics:mean"));
			statistics1_.setVariance(param_.getValue("statistics:variance1"));
			statistics2_.setVariance(param_.getValue("statistics:variance2"));
			
			setSamples();
		}

		void BiGaussModel::setOffset(double offset)
		{
			double diff = offset - getInterpolation().getOffset();
			min_ += diff;
			max_ += diff;
			statistics1_.setMean(statistics1_.mean()+diff);
			statistics2_.setMean(statistics2_.mean()+diff);
			
			InterpolationModel<>::setOffset(offset);

			param_.setValue("bounding_box:min", min_);
			param_.setValue("bounding_box:max", max_);
			param_.setValue("statistics:mean", statistics1_.mean());
		}

		const BiGaussModel::CoordinateType BiGaussModel::getCenter() const
		{
			return statistics2_.mean();
		}

}
