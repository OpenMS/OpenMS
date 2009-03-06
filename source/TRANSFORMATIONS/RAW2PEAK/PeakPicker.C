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
// $Maintainer: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPicker.h>

//Derived classes are included here
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>


namespace OpenMS
{
	//register products here
  void PeakPicker::registerChildren()
	{
		Factory<PeakPicker>::registerProduct ( PeakPickerCWT::getProductName(), &PeakPickerCWT::create );
		Factory<PeakPicker>::registerProduct ( PeakPickerHiRes::getProductName(), &PeakPickerHiRes::create );
	}

	PeakPicker::PeakPicker()
		: FactoryProduct("PeakPicker"),
			ProgressLogger()
	{
	}

	PeakPicker::~PeakPicker()
	{
	}
	
	void PeakPicker::pick(const MSSpectrum<>& /*input*/, MSSpectrum<>& /*output*/)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__); 
	}

	void PeakPicker::pickExperiment(const MSExperiment<>& /*input*/, MSExperiment<>& /*output*/)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__); 
	}

	void PeakPicker::updateMembers_()
	{
	}

} 
