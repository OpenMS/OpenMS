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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <numeric>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{

	IsotopeFitter1D::IsotopeFitter1D()
		: MaxLikeliFitter1D()
	{
		setName(getProductName());

		defaults_.setValue("statistics:variance",1.0,"Variance of the model.", StringList::create("advanced"));
		defaults_.setValue("charge",1,"Charge state of the model.", StringList::create("advanced"));
		defaults_.setValue("isotope:stdev",1.0,"Standard deviation of gaussian applied to the averagine isotopic pattern to simulate the inaccuracy of the mass spectrometer.", StringList::create("advanced"));
		defaults_.setValue("isotope:maximum",100,"Maximum isotopic rank to be considered.", StringList::create("advanced"));
		defaults_.setValue("interpolation_step",0.1,"Sampling rate for the interpolation of the model function.", StringList::create("advanced"));

		defaultsToParam_();
	}

	IsotopeFitter1D::IsotopeFitter1D(const IsotopeFitter1D& source)
		: MaxLikeliFitter1D(source)
	{
		updateMembers_();
	}

	IsotopeFitter1D::~IsotopeFitter1D()
	{
	}

	IsotopeFitter1D& IsotopeFitter1D::operator = (const IsotopeFitter1D& source)
	{
		if (&source == this) return *this;

		MaxLikeliFitter1D::operator = (source);
		updateMembers_();

		return *this;
	}

	IsotopeFitter1D::QualityType IsotopeFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
	{
		// Calculate bounding box
		min_ = max_ = set[0].getPos();
		for ( UInt pos=1; pos < set.size(); ++pos)
		{
			CoordinateType tmp = set[pos].getPos();
			if ( min_ > tmp ) min_ = tmp;
			if ( max_ < tmp ) max_ = tmp;
		}

		// Enlarge the bounding box by a few multiples of the standard deviation
		{
			stdev1_ = sqrt ( statistics_.variance() ) * tolerance_stdev_box_;
			min_ -= stdev1_;
			max_ += stdev1_;
		}

		// build model
		if (charge_==0)
		{
			model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("GaussModel"));
			model->setInterpolationStep( interpolation_step_ );

			Param tmp;
			tmp.setValue( "bounding_box:min", min_ );
			tmp.setValue( "bounding_box:max", max_ );
			tmp.setValue( "statistics:variance", statistics_.variance() );
			tmp.setValue( "statistics:mean", statistics_.mean() );
			model->setParameters( tmp );
		}
		else
		{
			model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("IsotopeModel"));

			Param iso_param = this->param_.copy( "isotope_model:", true );
			iso_param.removeAll("stdev");
			model->setParameters( iso_param );
			model->setInterpolationStep( interpolation_step_ );

			Param tmp;
			tmp.setValue( "statistics:mean", statistics_.mean() );
			tmp.setValue( "charge", static_cast<Int>( charge_ ) );
			tmp.setValue( "isotope:mode:GaussianSD", isotope_stdev_ );
			tmp.setValue( "isotope:maximum", max_isotope_ );

			model->setParameters( tmp );
      (static_cast<IsotopeModel*> (model)) -> setSamples( (static_cast<IsotopeModel*> (model))->getFormula() );
		}

		// fit offset
		QualityType quality;
		quality = fitOffset_(model, set, stdev1_, stdev1_, interpolation_step_);
		if (boost::math::isnan(quality) ) quality = -1.0;

		return quality;
	}

	void IsotopeFitter1D::updateMembers_()
	{
		MaxLikeliFitter1D::updateMembers_();
		statistics_.setVariance(param_.getValue("statistics:variance"));
		charge_ = param_.getValue("charge");
		isotope_stdev_ = param_.getValue("isotope:stdev");
		max_isotope_ = param_.getValue("isotope:maximum");
	}

}
