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
// $Maintainer: Andreas Bertsch $
// $Authors: $
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
    defaults_.setValue("window_size", 2.0, "The size of the m/z window where the peaks are removed, +/- window_size.");
		defaults_.setValue("default_charge", 2, "If the precursor has no charge set, the default charge is assumed.");
		defaults_.setValue("clean_all_charge_states", 1, "Set to 1 if precursor ions of all possible charge states should be removed.", StringList::create("advanced"));
		defaults_.setValue("consider_NH3_loss", 1, "Whether NH3 loss peaks from the precursor should be removed.");
		defaults_.setValue("consider_H2O_loss", 1, "Whether H2O loss peaks from the precursor should be removed.");
		defaults_.setValue("reduce_by_factor", 0, "Reduce the intensities of the precursor and related ions by a given factor (set 'set_to_zero' to 0).", StringList::create("advanced"));
		defaults_.setValue("factor", 1000.0, "Factor which is used to reduce the intensities if 'reduce_by_factor' is selected.", StringList::create("advanced"));
		defaults_.setValue("set_to_zero", 1, "Reduce the intensities of the precursor and related ions to zero.", StringList::create("advanced"));
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
