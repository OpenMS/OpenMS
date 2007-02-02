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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>

using namespace std;

namespace OpenMS
{

  ParentPeakMower::ParentPeakMower()
    : PreprocessingFunctor()
  {
		setName(ParentPeakMower::getProductName());
    defaults_.setValue("window_size", 2.0);
		defaults_.setValue("default_charge", (int)2);
		defaults_.setValue("clean_all_charge_states", (short)1);
		defaults_.setValue("consider_NH3_loss", (short)1);
		defaults_.setValue("consider_H2O_loss", (short)1);
		defaults_.setValue("reduce_by_factor", (short)0);
		defaults_.setValue("factor", 1000.0);
		defaults_.setValue("set_to_zero", (short)1);
		defaultsToParam_();
  }

  ParentPeakMower::ParentPeakMower(const ParentPeakMower& source)
    : PreprocessingFunctor(source)
  {
  }

  ParentPeakMower::~ParentPeakMower()
  {
  }

  ParentPeakMower& ParentPeakMower::operator = (const ParentPeakMower& source)
  {
		if (this != &source)
		{
    	PreprocessingFunctor::operator = (source);
		}
    return *this;
  }

	void ParentPeakMower::filterPeakSpectrum(PeakSpectrum& spectrum)
  {
    filterSpectrum(spectrum);
  }

  void ParentPeakMower::filterPeakMap(PeakMap& exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
