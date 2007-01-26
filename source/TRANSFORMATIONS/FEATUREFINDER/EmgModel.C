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



#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
#include <numeric>
#include <math.h>


namespace OpenMS
{
	EmgModel::EmgModel()
		: InterpolationModel()
	{
		setName(getProductName());
		
		defaults_.setValue("bounding_box:min",0.0f);
		defaults_.setValue("bounding_box:max",1.0f);
		defaults_.setValue("statistics:mean",0.0f);
		defaults_.setValue("statistics:variance",1.0f);
		defaults_.setValue("emg:height",100000.0f);
		defaults_.setValue("emg:width",5.0f);
		defaults_.setValue("emg:symmetry",5.0f);
		defaults_.setValue("emg:retention",1200.0f);
		
		defaultsToParam_();
	}

	EmgModel::EmgModel(const EmgModel& source)
		: InterpolationModel(source)
	{
		setParam(source.statistics_, source.height_, source.width_, source.symmetry_, source.retention_, source.min_, source.max_);
		updateMembers_();
	}

	EmgModel::~EmgModel()
	{	
	}

	EmgModel& EmgModel::operator = (const EmgModel& source)
	{
		if (&source == this) return *this;
		
		InterpolationModel::operator = (source);
		setParam(source.statistics_, source.height_, source.width_, source.symmetry_, source.retention_, source.min_, source.max_);
		updateMembers_();
		
		return *this;
	}

	void EmgModel::setSamples()
	{
		ContainerType& data = interpolation_.getData();
		data.clear();
		if (max_==min_) return;
		data.reserve( Size ( (max_-min_) / interpolation_step_ + 1 ) );
		CoordinateType pos = min_;

		double sqrt_2pi = sqrt(2*M_PI);
		double term_sq2 = (-2.4055/sqrt(2));
		double part1    = (height_*width_/symmetry_);
		double part2    = (width_*width_)/(2*symmetry_*symmetry_);
		double part3    = width_/symmetry_;

		for ( Size i = 0; pos< max_; ++i)
		{
			pos = min_ + i * interpolation_step_;

			// data.push_back (Simplified EMG)
			data.push_back((part1*sqrt_2pi*exp(part2-((pos-retention_)/symmetry_))/(1+exp(term_sq2*(((pos-retention_)/width_)-part3)))));
		}

		interpolation_.setScale  ( interpolation_step_ );
		interpolation_.setOffset ( min_ );
	}

	void EmgModel::setParam(const BasicStatistics& statistics, CoordinateType height, CoordinateType width, CoordinateType symmetry, CoordinateType retention, CoordinateType min, CoordinateType max)
	{
		statistics_.setMean(statistics.mean());
		statistics_.setVariance(statistics.variance());
		min_ = min;
		max_ = max;
		height_ = height;
		width_ =  width;
		symmetry_ =  symmetry;
		retention_ = retention;

		param_.setValue("bounding_box:min", min_);
		param_.setValue("bounding_box:max", max_);
		param_.setValue("statistics:mean", statistics_.mean());
		param_.setValue("statistics:variance", statistics_.variance());
		param_.setValue("emg:height", height_);
		param_.setValue("emg:width", width_);
		param_.setValue("emg:symmetry", symmetry_);
		param_.setValue("emg:retention", retention_);

		setSamples();
	}

	void EmgModel::setOffset(CoordinateType offset)
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

	const EmgModel::CoordinateType EmgModel::getCenter() const
	{
		return statistics_.mean();
	}
	
	void EmgModel::updateMembers_()
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
		
		setSamples();
	}
	
}
