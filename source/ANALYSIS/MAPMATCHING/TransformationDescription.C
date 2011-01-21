// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
  TransformationDescription::TransformationDescription():
		data_(TransformationDescription::DataPoints()), model_type_("none"), 
		model_(new TransformationModel())
  {
  }
	
  TransformationDescription::TransformationDescription(
		const TransformationDescription::DataPoints& data):
		data_(data), model_type_("none"), model_(new TransformationModel())
  {
  }

  TransformationDescription::~TransformationDescription()
  {
    delete model_;
  }

  TransformationDescription::TransformationDescription(
		const TransformationDescription& rhs)
  {
		data_ = rhs.data_;
		model_type_ = "none";
		model_ = 0; // initialize this before the "delete" call in "fitModel"!
		Param params;
		rhs.getModelParameters(params);
		fitModel(rhs.model_type_, params);
  }

  TransformationDescription& TransformationDescription::operator=(
		const TransformationDescription& rhs)
  {
    if (this == &rhs) return *this;

		data_ = rhs.data_;
		model_type_ = "none";
		Param params;
		rhs.getModelParameters(params);
		fitModel(rhs.model_type_, params);

    return *this;
  }

	void TransformationDescription::fitModel(const String& model_type, 
																					 const Param& params)
	{
		// if the transformation is the identity, don't fit another model:
		if (model_type_ == "identity") return;

		delete model_;
		model_ = 0; // avoid segmentation fault in case of exception
		if ((model_type == "none") || (model_type == "identity"))
		{
			model_ = new TransformationModel();
		}
		else if (model_type == "linear")
		{
			model_ = new TransformationModelLinear(data_, params);
			// // debug output:
			// DoubleReal slope, intercept;
			// TransformationModelLinear* lm = dynamic_cast<TransformationModelLinear*>(model_);
			// lm->getParameters(slope, intercept);
			// cout << "slope: " << slope << ", intercept: " << intercept << endl;
		}
		else if (model_type == "interpolated")
		{
			model_ = new TransformationModelInterpolated(data_, params);
		}
		else if (model_type == "b_spline")
		{
			model_ = new TransformationModelBSpline(data_, params);
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "unknown model type '" + model_type + "'");
		}
		model_type_ = model_type;
	}

	DoubleReal TransformationDescription::apply(DoubleReal value) const
	{
		return model_->evaluate(value);
	}

	const String& TransformationDescription::getModelType() const
	{
		return model_type_;
	}

	void TransformationDescription::getModelTypes(StringList& result)
	{
		result = StringList::create("linear,b_spline,interpolated");
		// "none" and "identity" don't count
	}

	void TransformationDescription::setDataPoints(const DataPoints& data)
	{
		data_ = data;
		model_type_ = "none"; // reset the model even if it was "identity"
		delete model_;
		model_ = new TransformationModel();
	}

	const TransformationDescription::DataPoints& 
	TransformationDescription::getDataPoints() const
	{
		return data_;
	}
	
  void TransformationDescription::getModelParameters(Param& params) const
	{
		model_->getParameters(params);
	}
		
	void TransformationDescription::invert()
	{
		for (TransformationDescription::DataPoints::iterator it = data_.begin();
				 it != data_.end(); ++it)
		{
			*it = make_pair(it->second, it->first);
		}
		// ugly hack for linear model with explicit slope/intercept parameters:
		if ((model_type_ == "linear") && data_.empty())
		{
			TransformationModelLinear *lm = 
				dynamic_cast<TransformationModelLinear*>(model_);
			lm->invert();
		}
		else
		{
			Param params;
			getModelParameters(params);
			fitModel(model_type_, params);
		}
	}

} // end of namespace OpenMS

