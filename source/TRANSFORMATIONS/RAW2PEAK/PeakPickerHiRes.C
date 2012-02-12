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
// $Maintainer: Erhan Kenar$
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <vector>

using namespace std;

namespace OpenMS
{
  PeakPickerHiRes::PeakPickerHiRes()
		: DefaultParamHandler("PeakPickerHiRes"),
			ProgressLogger()
  {
		// set default parameter values
                defaults_.setValue("signal_to_noise", 1.0, "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)");
		defaults_.setMinFloat("signal_to_noise", 0.0);

		defaults_.setValue("ms1_only", "false", "If true, peak picking is only applied to MS1 scans. Other scans are copied to the output without changes.");
		defaults_.setValidStrings("ms1_only", StringList::create("true,false"));
		
		// parameters for SNTestimator config ...

		// write defaults into Param object param_
		defaultsToParam_();

		
		// initialize class members
		signal_to_noise_ = param_.getValue("signal_to_noise");
  }

  PeakPickerHiRes::~PeakPickerHiRes()
  {
  }

	void PeakPickerHiRes::updateMembers_()
	{
		signal_to_noise_ = param_.getValue("signal_to_noise");
	}
}
