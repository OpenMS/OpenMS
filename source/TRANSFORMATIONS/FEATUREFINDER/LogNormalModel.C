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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>
#include <numeric>
#include <math.h>

namespace OpenMS
{
	LogNormalModel::LogNormalModel()
		: InterpolationModel()
	{
		setName(getProductName());
		
		defaults_.setValue("bounding_box:min",0.0f,"lower bound of bounding box");
		defaults_.setValue("bounding_box:max",1.0f,"upper bound of bounding box");
		defaults_.setValue("statistics:mean",0.0f,"mean");
		defaults_.setValue("statistics:variance",1.0f,"variance");
		defaults_.setValue("emg:height",100000.0f,"height");
		defaults_.setValue("emg:width",5.0f,"width");
		defaults_.setValue("emg:symmetry",5.0f,"symmetry factor");
		defaults_.setValue("emg:retention",1200.0f,"retention");
		defaults_.setValue("lognormal:r",2.0f,"lognormal scale");

		defaultsToParam_();
	}

	LogNormalModel::LogNormalModel(const LogNormalModel& source)
		: InterpolationModel(source)
	{
		setParameters( source.getParameters() );
		updateMembers_();
	}

	LogNormalModel::~LogNormalModel()
	{
	}

	LogNormalModel& LogNormalModel::operator = (const LogNormalModel& source)
	{
		if (&source == this) return *this;
		
		InterpolationModel::operator = (source);
		setParameters( source.getParameters() );
		updateMembers_();
		
		return *this;
	}

	void LogNormalModel::setSamples()
	{
		ContainerType& data = interpolation_.getData();
		data.clear();
		if (max_==min_) return;
		data.reserve( UInt ( (max_-min_) / interpolation_step_ + 1 ) );
		CoordinateType pos = min_;

		double canelValue = retention_ - (width_*symmetry_)/(symmetry_*symmetry_ - 1);
		for ( UInt i = 0; pos< max_; ++i)
		{
			pos = min_ + i * interpolation_step_;

			if (pos <= canelValue)
			data.push_back(0);
			else
				data.push_back( height_*exp(-log(r_)/(log(symmetry_)*log(symmetry_))*log((pos-retention_)*(symmetry_*symmetry_-1)/width_/symmetry_+1)*log((pos-retention_)*(symmetry_*symmetry_-1)/width_/symmetry_+1)));
		}

		interpolation_.setScale  ( interpolation_step_ );
		interpolation_.setOffset ( min_ );
	}

	void LogNormalModel::setOffset(CoordinateType offset)
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

	LogNormalModel::CoordinateType LogNormalModel::getCenter() const
	{
		return statistics_.mean();
	}

	void  LogNormalModel::updateMembers_()
	{
		InterpolationModel::updateMembers_();
		
		min_ = param_.getValue("bounding_box:min");
		max_ = param_.getValue("bounding_box:max");
		statistics_.setMean( param_.getValue("statistics:mean") );
		statistics_.setVariance(param_.getValue("statistics:variance"));
		height_ = param_.getValue("emg:height");
		width_ = param_.getValue("emg:width");
		symmetry_ = param_.getValue("emg:symmetry");
		retention_ = param_.getValue("emg:retention");
		r_ = param_.getValue("lognormal:r");
		
		setSamples();
	}

}
