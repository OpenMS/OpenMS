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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>

using namespace std;
namespace OpenMS
{

  NLargest::NLargest()
    : PreprocessingFunctor()
  {
		setName(NLargest::getProductName());
    defaults_.setValue("n", 200, "The number of peaks to keep.");
		defaultsToParam_();
  }

	NLargest::NLargest(UInt n)
		: PreprocessingFunctor()
	{
		setName(NLargest::getProductName());
		defaults_.setValue("n", n, "The number of peaks to keep");
		defaultsToParam_();
	}

  NLargest::NLargest(const NLargest& source)
    : PreprocessingFunctor(source)
  {
  }

  NLargest::~NLargest()
  {
  }

  NLargest& NLargest::operator=(const NLargest& source)
  {
		if (this != &source)
		{
    	PreprocessingFunctor::operator=(source);
		}
    return *this;
  }

  void NLargest::filterPeakSpectrum(PeakSpectrum& spectrum)
  {
    filterSpectrum(spectrum);
  }

  void NLargest::filterPeakMap(PeakMap& exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }


}
