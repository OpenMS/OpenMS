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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/FeatureEditDialog.h>

using namespace std;

namespace OpenMS
{

	FeatureEditDialog::FeatureEditDialog(QWidget* parent)
		: QDialog(parent),
			feature_()
	{
		setupUi(this);
	}

	void FeatureEditDialog::setFeature(const Feature& feature)
	{
		//copy feature
		feature_ = feature;
		//update widgets according to feature data
		mz_->setValue(feature_.getMZ());
		rt_->setValue(feature_.getRT());
		int_->setValue(feature_.getIntensity());
		charge_->setValue(feature_.getCharge());
	}
	
	const Feature& FeatureEditDialog::getFeature() const
	{
		//update feature data according to widget
		feature_.setMZ(mz_->value());
		feature_.setRT(rt_->value());
		feature_.setIntensity(int_->value());
		feature_.setCharge(charge_->value());
		
		//return feature
		return feature_;
	}
	
} // namespace
