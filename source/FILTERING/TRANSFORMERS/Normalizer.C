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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <cmath>
#include <set>

using namespace std;
namespace OpenMS
{
	Normalizer::Normalizer()
		:DefaultParamHandler("Normalizer")
  {	
    defaults_.setValue("method", "to_one", "Normalize by deviding though the TIC ('to_TIC') or normalize to max intensity of one ('to_one').");
    defaults_.setValidStrings("method", StringList::create("to_one,to_TIC"));
		defaultsToParam_();
  }
	
	Normalizer::~Normalizer()
  {
  }

  Normalizer::Normalizer(const Normalizer& source)
    : DefaultParamHandler(source)
  {
  }

  Normalizer& Normalizer::operator=(const Normalizer& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }
	
  void Normalizer::filterPeakSpectrum(PeakSpectrum& spectrum)
  {
    filterSpectrum(spectrum);
  }

  void Normalizer::filterPeakMap(PeakMap& exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }
	
}
