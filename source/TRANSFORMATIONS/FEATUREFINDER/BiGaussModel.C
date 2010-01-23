// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>
#include <numeric>

namespace OpenMS
{
    BiGaussModel::BiGaussModel()
		: InterpolationModel(), statistics1_(), statistics2_()
		{
			setName(getProductName());

      defaults_.setValue("bounding_box:min",0.0,"Lower end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("bounding_box:max",1.0,"Upper end of bounding box enclosing the data used to fit the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:mean",0.0,"Centroid position of the model, this also separates both halves of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance1",1.0,"Variance of the first gaussian, used for the lower half of the model.", StringList::create("advanced"));
      defaults_.setValue("statistics:variance2",1.0,"Variance of the second gaussian, used for the upper half of the model.", StringList::create("advanced"));

			defaultsToParam_();
		}

  	BiGaussModel::BiGaussModel(const BiGaussModel& source)
		: InterpolationModel(source)
		{
			setParameters( source.getParameters() );
			updateMembers_();
		}

    BiGaussModel::~BiGaussModel()
    {
    }

   	BiGaussModel& BiGaussModel::operator = (const BiGaussModel& source)
		{
			if (&source == this) return *this;

			InterpolationModel::operator = (source);
			setParameters( source.getParameters() );
			updateMembers_();

			return *this;
		}

		void BiGaussModel::setSamples()
		{
			LinearInterpolation::container_type& data = interpolation_.getData();
			data.clear();
			if (max_==min_) return;
			data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
			CoordinateType pos = min_;

			for ( UInt i = 0; pos< max_; ++i)
			{
				pos = min_ + i * interpolation_step_;
				if (pos < statistics1_.mean())
        {
					data.push_back( statistics1_.normalDensity_sqrt2pi(pos) );
        }
				else
        {
					data.push_back( statistics2_.normalDensity_sqrt2pi(pos) );
        }
			}
			// scale data so that integral over distribution equals one
			// multiply sum by interpolation_step_ -> rectangular approximation of integral
			IntensityType factor = scaling_ / interpolation_step_ /
						std::accumulate ( data.begin(), data.end(), IntensityType(0) );

			for (LinearInterpolation::container_type::iterator it=data.begin();	it!=data.end();	++it)
      {
				*it *= factor;
			}

			interpolation_.setScale  ( interpolation_step_ );
			interpolation_.setOffset ( min_ );
		}

		void BiGaussModel::updateMembers_()
		{
			InterpolationModel::updateMembers_();

			min_ = param_.getValue("bounding_box:min");
			max_ = param_.getValue("bounding_box:max");
			statistics1_.setMean(param_.getValue("statistics:mean"));
			statistics2_.setMean(param_.getValue("statistics:mean"));
			statistics1_.setVariance(param_.getValue("statistics:variance1"));
			statistics2_.setVariance(param_.getValue("statistics:variance2"));

			setSamples();
		}

    void BiGaussModel::setOffset(CoordinateType offset)
		{
			DoubleReal diff = offset - getInterpolation().getOffset();
			min_ += diff;
			max_ += diff;
			statistics1_.setMean(statistics1_.mean()+diff);
			statistics2_.setMean(statistics2_.mean()+diff);

			InterpolationModel::setOffset(offset);

			param_.setValue("bounding_box:min", min_);
			param_.setValue("bounding_box:max", max_);
			param_.setValue("statistics:mean", statistics1_.mean());
		}

		BiGaussModel::CoordinateType BiGaussModel::getCenter() const
		{
			return statistics2_.mean();
		}

}
