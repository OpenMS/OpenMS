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
// $Maintainer: Chris Bielow, Andreas Bertsch $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>

using namespace std;
namespace OpenMS
{

  SpectraMerger::SpectraMerger()
    : DefaultParamHandler("SpectraMerger")
  {
		defaults_.setValue("ms_levels", IntList::create("1"), "Spectra of that MS levels are subject to be merged.");
		defaults_.setMinInt("ms_levels", 1);

    // common
    defaults_.setValue("mz_binning_width", 10e-5, "Max m/z distance of two peaks to merged.", StringList::create("advanced"));
    defaults_.setMinFloat("mz_binning_width", 0);

    defaults_.setValue("mz_binning_width_units", "Da", "Unit in which the distance between two peaks is given.", StringList::create("advanced"));
    defaults_.setValidStrings("mz_binning_width_units", StringList::create("Da,ppm"));

    // block merging
    defaults_.setValue("block_method:rt_block_size", 5, "Number of scans to be summed up.");
    defaults_.setMinInt("block_method:rt_block_size", 1);

    // same precursor MS/MS merging
   	defaults_.setValue("precursor_method:mz_tolerance", 10e-5, "Max m/z distance of the precursor entries of two spectra to be merged in Dalton.");
		defaults_.setMinFloat("precursor_method:mz_tolerance", 0);
    defaults_.setValue("precursor_method:rt_tolerance", 5.0, "Max RT distance of the precursor entries of two spectra to be merged in Dalton.");
		defaults_.setMinFloat("precursor_method:rt_tolerance", 0);

		defaultsToParam_();
  }

  SpectraMerger::SpectraMerger(const SpectraMerger& source)
    : DefaultParamHandler(source)
  {
  }

  SpectraMerger::~SpectraMerger()
  {
  }

  SpectraMerger& SpectraMerger::operator=(const SpectraMerger& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }


}
