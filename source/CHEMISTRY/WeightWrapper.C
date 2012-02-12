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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/WeightWrapper.h>

namespace OpenMS{
	
	WeightWrapper::WeightWrapper()
	:weight_mode_(MONO)
	{
	}

	WeightWrapper::WeightWrapper(const WEIGHTMODE weight_mode)
	:weight_mode_(weight_mode)
	{
	}

	WeightWrapper::WeightWrapper(const WeightWrapper& source)
	 : weight_mode_(source.weight_mode_) 
	{
	}
	
	WeightWrapper::~WeightWrapper()
	{
	}
	
	void WeightWrapper::setWeightMode(const WEIGHTMODE mode)
	{
		if (mode >= WeightWrapper::SIZE_OF_WEIGHTMODE) throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "setWeightMode() received illegal 'mode' value!");
		weight_mode_ = mode;
	}

	WeightWrapper::WEIGHTMODE WeightWrapper::getWeightMode() const
	{
		return weight_mode_;
	}

	DoubleReal WeightWrapper::getWeight(const AASequence& aa) const
	{
		if (weight_mode_==WeightWrapper::MONO) return aa.getMonoWeight();
		else return aa.getAverageWeight();
	}
	
	DoubleReal WeightWrapper::getWeight(const EmpiricalFormula& ef) const
	{
		if (weight_mode_==WeightWrapper::MONO) return ef.getMonoWeight();
		else return ef.getAverageWeight();
	}
	

	DoubleReal WeightWrapper::getWeight(const Residue& r, Residue::ResidueType res_type) const
	{
		if (weight_mode_==WeightWrapper::MONO) return r.getMonoWeight(res_type);
		else return r.getAverageWeight(res_type);
	}

}
