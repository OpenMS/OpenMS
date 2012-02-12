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
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>

using namespace std;

namespace OpenMS
{

  NLargest::NLargest()
    : DefaultParamHandler("NLargest")
  {
    init_();
  }

  NLargest::NLargest(UInt n)
    : DefaultParamHandler("NLargest")
  {
    init_();
    // after initialising with the default value, use the provided n
    param_.setValue("n",n);
    updateMembers_();
  }

  void NLargest::init_()
  {
    defaults_.setValue("n", 200, "The number of peaks to keep");
    defaultsToParam_();
  }

  NLargest::~NLargest()
  {
  }

  NLargest::NLargest(const NLargest& source)
    : DefaultParamHandler(source)
  {
    updateMembers_();
  }

  NLargest& NLargest::operator=(const NLargest& source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
      updateMembers_();
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

  void NLargest::updateMembers_()
  {
    peakcount_ = (UInt)param_.getValue("n");
  }
}
